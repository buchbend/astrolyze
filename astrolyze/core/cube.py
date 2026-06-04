"""The PPV wrapper: :class:`Cube`, a thin composition over :class:`spectral_cube.SpectralCube`.

The heavy lifting — moments, the spectral axis, slicing — is delegated to ``spectral-cube``;
``Cube`` adds the astrolyze context (beam + rest frequency + velocity convention, sourced from
:class:`~astrolyze.io.Metadata`) and the type transitions that carry it: ``moment0() -> Map``
and ``cube[:, y, x] -> Spectrum`` (ADR-0004). It is **not** a ``SpectralCube`` subclass — it
holds one (``._sc``) and delegates, so we stay loosely coupled to upstream internals (ADR-0003).
"""

from __future__ import annotations

import numpy as np
import astropy.units as u
import radio_beam
from astropy.convolution import Gaussian1DKernel
from radio_beam.utils import BeamError
from spectral_cube import SpectralCube

from astrolyze.io import Metadata
from astrolyze.units import (
    BrightnessTemperatureScale,
    CalibrationScale,
    MissingContextError,
    coerce_calibration_scale,
    coerce_temperature_scale,
    convert,
)

from . import _coords
from ._base import ContextCarrier, _emit
from .map import Map
from .spectrum import Spectrum

#: The authoritative surface-brightness representation a harmonised cube lands in (I_ν).
_SURFACE_BRIGHTNESS_UNIT = u.MJy / u.sr

#: Gaussian FWHM <-> sigma. Spectral smoothing kernels are specified by sigma; users state a
#: resolution (a width = FWHM), which is the natural, observable quantity.
_FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


class LossyDirectionError(ValueError):
    """Raised when a matching op is asked to *gain* resolution rather than degrade it.

    The matching ops (:meth:`Cube.convolve_to_beam`, :meth:`Cube.spectral_bin`,
    :meth:`Cube.spectral_smooth_to`) only ever *lose* resolution — they smooth to a larger
    beam or coarser channels. The inverse (deconvolving to a smaller beam, up-sampling to
    finer channels) would invent structure the data never contained; astrolyze refuses it
    rather than super-resolving silently (no silent physics, ADR-0003). A ``ValueError`` so
    it reads as the programming error it is.
    """


class Cube(ContextCarrier):
    """A PPV cube: a wrapped :class:`~spectral_cube.SpectralCube` that carries context."""

    _viz_function = "plot_cube"

    def __init__(self, spectral_cube: SpectralCube, metadata: Metadata):
        self._sc = spectral_cube
        self.metadata = metadata

    @classmethod
    def from_loaded(cls, loaded) -> "Cube":
        """Build a :class:`Cube` from an :class:`~astrolyze.io.LoadedData` (the io seam).

        The beam is attached at construction rather than via ``with_beam`` — spectral-cube
        refuses ``with_beam`` on a Jy/beam cube (it would change the spatial resolution), but
        we are only recording the beam the header already states, not smoothing.
        """
        bunit = loaded.metadata.bunit or u.dimensionless_unscaled
        # A dask-backed source (the lazy Zarr backend, #23) must stay a dask array: wrapping it
        # in u.Quantity materialises it to numpy, so we keep it bare and carry the unit via the
        # cube meta instead (see _build_cube). FITS data is a plain ndarray -> the eager path.
        data = loaded.data if _is_dask(loaded.data) else u.Quantity(loaded.data, bunit)
        sc = _build_cube(data, loaded.wcs, loaded.metadata.beam, unit=bunit)
        return cls(sc, loaded.metadata)

    # -- Zarr backend (issue #23): same context as the FITS path, via the io seam -------
    @classmethod
    def from_zarr(cls, store) -> "Cube":
        """Build a :class:`Cube` from a Zarr store, carrying the same context as the FITS path.

        A thin convenience over the io seam: ``Cube.from_loaded(load(store))``. ``load``
        dispatches the Zarr store to the lazy backend, and :meth:`from_loaded` attaches the
        beam / rest frequency / convention exactly as it does for FITS (ADR-0004/0006)."""
        from astrolyze.io import load

        return cls.from_loaded(load(store))

    def to_zarr(self, directory, **layout):
        """Write this cube to an xarray-native Zarr v3 store under *directory*; return its path.

        A thin convenience over the io seam (``io.save(..., format="zarr")``): the schema, the
        verbatim FITS-WCS string, and the data are carried unchanged, so a later
        :meth:`from_zarr` reconstructs the same context. ``**layout`` (``chunks`` / ``shards`` /
        ``compressors``) is the caller's Zarr layout, passed straight through — astrolyze fixes
        no chunking policy."""
        from astrolyze.io import save

        return save(
            self._data_quantity.value,
            self.metadata,
            directory,
            format="zarr",
            base_header=self._sc.wcs.to_header(),
            **layout,
        )

    # -- delegated PPV operations (type transitions carry context) ----------------------
    def moment0(self) -> Map:
        """The zeroth moment (velocity-integrated intensity) as a :class:`Map`."""
        return self.moment(order=0)

    def moment(self, order: int = 0, axis: int = 0) -> Map:
        """Delegate the moment to spectral-cube and wrap the 2D result, carrying context.

        The moment changes the physical unit (e.g. Jy/beam -> Jy/beam km/s); the beam, rest
        frequency and convention are preserved on the resulting :class:`Map`."""
        proj = self._sc.moment(order=order, axis=axis)
        _emit("moment", params={"order": order, "axis": axis})
        return Map(
            proj, getattr(proj, "wcs", None), self._metadata_with_unit(proj.unit)
        )

    def __getitem__(self, key):
        """Slice the underlying cube and wrap by resulting dimensionality, carrying context:
        3D -> :class:`Cube`, 2D channel -> :class:`Map`, 1D spectrum -> :class:`Spectrum`.
        """
        result = self._sc[key]
        if isinstance(result, SpectralCube):
            return Cube(result, self.metadata)
        ndim = getattr(result, "ndim", None)
        if ndim == 1:  # spectral_cube.OneDSpectrum
            return Spectrum.from_oned(result, self.metadata)
        if ndim == 2:  # a channel map (slicing does not change the unit)
            return Map(
                result,
                getattr(result, "wcs", None),
                self._metadata_with_unit(result.unit),
            )
        return result  # 0D / anything else: hand back the raw upstream value

    # -- guarded beam / channel matching (issue #31) ------------------------------------
    # astrolyze stays thin here: spectral-cube does the convolution/binning maths and
    # radio_beam decides which direction is resolution-losing. The only value added is the
    # lossy-direction guard (never super-resolve, ADR-0003) and the context carry (the new,
    # larger beam recorded on the metadata, ADR-0004). Noise propagation is issue #32.
    def convolve_to_beam(self, beam: radio_beam.Beam) -> "Cube":
        """Smooth the cube spatially to a **larger** *beam* via spectral-cube ``convolve_to``.

        The target must be larger than the current beam in the sense that it can be reached
        by convolution — i.e. ``beam`` deconvolved by the current beam exists. A smaller,
        equal, or otherwise non-degrading target would require deconvolution / super-
        resolution and so :class:`LossyDirectionError` is raised instead (ADR-0003).

        Returns a new :class:`Cube` carrying the new, larger beam as context (ADR-0004)."""
        self._require_larger_beam(beam)
        smoothed = self._masked_sc().convolve_to(beam)
        _emit("convolve_to_beam", params={"beam": str(beam)})
        return Cube(smoothed, self._metadata_with_beam(beam))

    def spectral_bin(self, factor: int) -> "Cube":
        """Bin the spectral axis by an integer *factor* via spectral-cube ``downsample_axis``.

        Coarsens the spectral resolution (fewer, wider channels). A *factor* of 1 (no-op) or
        less (which would imply finer channels / up-sampling) raises
        :class:`LossyDirectionError` — ``spectral_bin`` only ever coarsens (ADR-0003).

        Returns a new :class:`Cube`; the spatial beam is unchanged and carried through."""
        if not isinstance(factor, (int, np.integer)) or factor <= 1:
            raise LossyDirectionError(
                f"spectral_bin only coarsens: factor must be an integer > 1, got {factor!r} "
                "(a finer/up-sampled grid is never invented)"
            )
        binned = self._masked_sc().downsample_axis(int(factor), axis=0)
        _emit("spectral_bin", params={"factor": int(factor)})
        return Cube(binned, self._metadata_with_beam(self.metadata.beam))

    def spectral_smooth_to(self, width: u.Quantity) -> "Cube":
        """Smooth the spectral axis to a **broader** resolution *width* (a FWHM).

        Delegates the Gaussian smoothing to spectral-cube ``spectral_smooth``. The target
        *width* must be broader than the native channel width; a finer target would be up-
        sampling and raises :class:`LossyDirectionError` (ADR-0003). The convolving kernel is
        the one whose FWHM, added in quadrature to the channel width, reaches *width*.

        Returns a new :class:`Cube`; the channel grid and the spatial beam are unchanged and
        carried through."""
        channel_width = self._channel_width()
        target = width.to(channel_width.unit, equivalencies=u.spectral())
        if target <= channel_width:
            raise LossyDirectionError(
                f"spectral_smooth_to only degrades: target width {width} is not broader than "
                f"the native channel width {channel_width} (astrolyze never up-samples to a "
                "finer spectral resolution)"
            )
        # Quadrature: the kernel that takes the native channel width up to the target FWHM.
        kernel_fwhm = np.sqrt(target**2 - channel_width**2)
        sigma_channels = (kernel_fwhm / channel_width).to_value(
            u.dimensionless_unscaled
        ) / _FWHM_PER_SIGMA
        smoothed = self._masked_sc().spectral_smooth(Gaussian1DKernel(sigma_channels))
        _emit("spectral_smooth_to", params={"width": str(width)})
        return Cube(smoothed, self._metadata_with_beam(self.metadata.beam))

    def match_to(
        self, other: "Cube", *, reproject: bool = False
    ) -> tuple["Cube", "Cube"]:
        """Bring ``self`` and ``other`` to a common beam (a common-beam + line-ratio helper).

        The common beam is the smallest beam both cubes can be *smoothed* to (radio_beam's
        common beam) — never smaller than either, so neither cube is super-resolved. Each
        cube is convolved to it (the one already at the common beam is left untouched).

        Reproject is **explicit, never a side effect**: only when ``reproject=True`` is
        ``self`` reprojected onto ``other``'s spatial grid so the pair shares one pixel grid
        (their spectral axes are left alone). Line-ratio work must opt in to regridding.

        Returns ``(matched_self, matched_other)`` as :class:`Cube`\\ s carrying the new
        common beam (ADR-0004)."""
        common = radio_beam.commonbeam.commonbeam(
            radio_beam.Beams(beams=[self.metadata.beam, other.metadata.beam])
        )
        matched_self = self._convolve_if_needed(common)
        matched_other = other._convolve_if_needed(common)
        if reproject:
            matched_self = matched_self._reproject_spatial_to(matched_other)
        _emit("match_to", params={"beam": str(common), "reproject": reproject})
        return matched_self, matched_other

    # -- surface-brightness harmonisation (issue #25) -----------------------------------
    # I_ν (MJy/sr) is the one representation every input maps onto, because specific
    # intensity is beam- and calibration-scale-independent (a physical quantity of the sky,
    # not of the instrument). The new physics here is the calibration TEMPERATURE SCALE: a
    # Kelvin cube must DECLARE whether it is T_mb / T_A* / T_R* before we may treat it as a
    # brightness temperature (K is not silently T_mb, ADR-0003), and T_A* additionally needs
    # eta_mb to reach the main-beam scale. The unit maths is reused from units.convert; this
    # adds only that calibration gate plus the per-channel evaluation.
    def harmonize_to_surface_brightness(self) -> "Cube":
        """Harmonise this cube to specific intensity I_ν in MJy/sr — the authoritative,
        instrument-independent representation (issue #25, ADR-0003/0004).

        The conversion is routed by the source unit:

        - **temperature (K)** — requires a declared :attr:`calibration_scale` (``T_mb`` /
          ``T_A*`` / ``T_R*``); ``K`` is never silently assumed to be ``T_mb``. ``T_A*`` also
          requires :attr:`eta_mb`, applied as ``T_mb = T_A*/eta_mb`` before the brightness
          conversion. K↔I_ν is the **per-channel** Rayleigh-Jeans relation
          ``I_ν = 2kν²/c²·T_mb`` evaluated at each channel's own absolute frequency (not the
          band centre);
        - **per-beam flux (Jy/beam)** — pure beam geometry via the beam solid angle;
        - **surface brightness (Jy/sr / MJy/sr)** — a plain unit rescale.

        Returns a new :class:`Cube` in MJy/sr carrying this object's context (ADR-0004); the
        calibration scale + efficiency are recorded so the derived K view is reversible. Raises
        :class:`~astrolyze.units.MissingContextError` for any absent required context."""
        src = self.unit
        if _is_temperature_unit(src):
            data = self._main_beam_temperature()  # applies T_A*/eta_mb when needed
            # K -> I_ν is the per-channel RJ relation: I_ν = factor(ν)·T_mb. Multiplying by
            # the per-channel factor (and dividing by it in as_kelvin) is the exact inverse
            # pair, so the K -> I_ν -> K round trip is bit-stable per channel (issue #25).
            factor = self._per_channel_rj_factor()
            harmonized = (data.to_value(u.K) * factor) * _SURFACE_BRIGHTNESS_UNIT
        else:
            # Jy/beam (geometry) or Jy/sr->MJy/sr (rescale): frequency-independent, so the
            # whole array converts in one call with the object's beam supplied where needed.
            harmonized = convert(
                self._data_quantity,
                _SURFACE_BRIGHTNESS_UNIT,
                beam=self.metadata.beam,
            )
        _emit(
            "harmonize_to_surface_brightness",
            params={"unit": str(_SURFACE_BRIGHTNESS_UNIT)},
        )
        return Cube(
            _build_cube(harmonized, self._sc.wcs, self.metadata.beam),
            self._metadata_with_unit(_SURFACE_BRIGHTNESS_UNIT),
        )

    def as_kelvin(self, *, temperature_scale=None) -> "Cube":
        """Derived **view** of the authoritative I_ν cube as main-beam brightness temperature.

        A non-mutating accessor (issue #25): it converts I_ν -> K per channel at each channel's
        own absolute frequency and returns a *new* :class:`Cube`, leaving this I_ν cube
        untouched. ``temperature_scale`` (the RJ-vs-Planck brightness *law*) is mandatory and
        never defaulted — it is the genuinely-ambiguous choice the cube cannot make for you
        (ADR-0003); omitting it raises :class:`~astrolyze.units.MissingContextError`."""
        if temperature_scale is None:
            raise MissingContextError(
                "temperature_scale is required (rayleigh_jeans | planck) for the Kelvin view: "
                "RJ-vs-Planck is the top silent-error trap in radio/sub-mm work, so astrolyze "
                "never assumes it"
            )
        intensity = self._data_quantity.to(_SURFACE_BRIGHTNESS_UNIT)
        if (
            coerce_temperature_scale(temperature_scale)
            is BrightnessTemperatureScale.RAYLEIGH_JEANS
        ):
            # RJ is linear: divide by the per-channel factor used on the way in — the exact
            # inverse of harmonisation, so I_ν -> K -> I_ν is bit-stable per channel.
            factor = self._per_channel_rj_factor()
            data = (intensity.to_value(_SURFACE_BRIGHTNESS_UNIT) / factor) * u.K
        else:
            # Planck is nonlinear; convert per channel at its own absolute frequency.
            data = self._per_channel_convert(intensity, u.K, temperature_scale="planck")
        return Cube(
            _build_cube(data, self._sc.wcs, self.metadata.beam),
            self._metadata_with_unit(u.K),
        )

    def as_jy_per_beam(self) -> "Cube":
        """Derived **view** of the authoritative I_ν cube as per-beam flux (Jy/beam).

        A non-mutating accessor (issue #25): pure beam geometry (I_ν is frequency-independent
        here), returning a *new* :class:`Cube` and leaving this I_ν cube untouched. Requires a
        beam; without one it raises rather than guess a beam shape (ADR-0003)."""
        data = convert(self._data_quantity, u.Jy / u.beam, beam=self.metadata.beam)
        return Cube(
            _build_cube(data, self._sc.wcs, self.metadata.beam),
            self._metadata_with_unit(u.Jy / u.beam),
        )

    # -- harmonisation helpers ----------------------------------------------------------
    def _main_beam_temperature(self) -> u.Quantity:
        """The data as main-beam brightness temperature T_mb, demanding the calibration scale.

        K is never silently T_mb: the cube must declare :attr:`calibration_scale`. ``T_A*`` is
        divided by :attr:`eta_mb` (required) to reach the main-beam scale; ``T_mb`` / ``T_R*``
        are already on it (ADR-0003)."""
        scale = self.metadata.calibration_scale
        if scale is None:
            raise MissingContextError(
                "a calibration_scale (T_mb | T_A* | T_R*) is required to treat this Kelvin "
                "cube as a brightness temperature; astrolyze never silently assumes K = T_mb "
                "(ADR-0003) — declare the scale on the header/metadata"
            )
        scale = coerce_calibration_scale(scale)
        data = self._data_quantity
        if scale is CalibrationScale.T_A_STAR:
            if self.metadata.eta_mb is None:
                raise MissingContextError(
                    "eta_mb (main-beam efficiency) is required for a T_A* cube: T_mb = "
                    "T_A*/eta_mb must be applied before converting to surface brightness; "
                    "astrolyze never assumes an efficiency"
                )
            data = data / float(self.metadata.eta_mb)
        return data

    def _per_channel_rj_factor(self) -> np.ndarray:
        """The per-channel Rayleigh-Jeans factor ``I_ν/T_mb`` (MJy/sr per K), one per channel.

        K↔I_ν depends on ν, and a PPV cube spans a range of frequencies, so a single
        band-centre frequency would bias every off-centre channel. Each factor is evaluated at
        that channel's **own** authoritative absolute frequency (``coordinates.frequency``,
        derived under the object's stated convention + rest frequency — absent, it raises;
        ADR-0003). It is shaped to broadcast over the (y, x) plane, so ``T·factor`` (harmonise)
        and ``I_ν/factor`` (the K view) are an exact inverse pair — bit-stable per channel."""
        frequencies = (
            self.coordinates.frequency
        )  # authoritative; raises if context absent
        factor = np.array(
            [
                convert(
                    1.0 * u.K,
                    _SURFACE_BRIGHTNESS_UNIT,
                    rest_frequency=nu,
                    temperature_scale="rayleigh_jeans",
                ).to_value(_SURFACE_BRIGHTNESS_UNIT)
                for nu in frequencies
            ],
            dtype="float64",
        )
        return factor[:, np.newaxis, np.newaxis]  # broadcast over (channel, y, x)

    def _per_channel_convert(self, data, target, *, temperature_scale):
        """Convert ``data`` (a 3D Quantity) to ``target`` channel by channel at each channel's
        own absolute frequency — the general (e.g. nonlinear Planck) path. Delegates the
        brightness maths to ``units.convert``; only the per-channel ν is added here."""
        frequencies = self.coordinates.frequency
        out = np.empty(data.shape, dtype="float64")
        target_unit = u.Unit(target)
        for k, nu in enumerate(frequencies):
            out[k] = convert(
                data[k],
                target_unit,
                rest_frequency=nu,
                temperature_scale=temperature_scale,
            ).to_value(target_unit)
        return out * target_unit

    # -- matching helpers ---------------------------------------------------------------
    def _require_larger_beam(self, beam: radio_beam.Beam) -> None:
        """Raise unless *beam* is reachable from the current beam by convolution (i.e. larger).

        radio_beam's ``deconvolve`` succeeds exactly when ``beam`` can be written as the
        current beam convolved with a real kernel — the test for "strictly larger". An equal
        or smaller target fails to deconvolve and is on the refused side of the guard."""
        try:
            beam.deconvolve(self.metadata.beam)
        except BeamError as exc:
            raise LossyDirectionError(
                f"target beam {beam} is not larger than the current beam "
                f"{self.metadata.beam}: convolving to it would require deconvolution / "
                "super-resolution, which astrolyze refuses (never invents structure)"
            ) from exc

    def _convolve_if_needed(self, beam: radio_beam.Beam) -> "Cube":
        """Convolve to *beam* unless already at it (used by :meth:`match_to`).

        Unlike :meth:`convolve_to_beam`, reaching the common beam may be a no-op for the
        coarser cube (it is already there), which must not raise."""
        try:
            self._require_larger_beam(beam)
        except LossyDirectionError:
            return (
                self  # already at (or effectively at) the common beam — nothing to do.
            )
        return self.convolve_to_beam(beam)

    def _reproject_spatial_to(self, other: "Cube") -> "Cube":
        """Reproject onto *other*'s spatial grid (a common grid), keeping our spectral axis."""
        reprojected = self._masked_sc().reproject(other._sc.header)
        return Cube(reprojected, self._metadata_with_beam(self.metadata.beam))

    def _channel_width(self) -> u.Quantity:
        """The native spectral channel width (a positive Quantity in the axis' units)."""
        axis = self._sc.spectral_axis
        return abs(axis[1] - axis[0])

    def _masked_sc(self) -> SpectralCube:
        """The underlying cube with a finite mask attached.

        spectral-cube's spatial operations (``convolve_to``, ``reproject``) require a mask;
        a cube built straight from an array has none. We attach an all-finite mask so the
        delegation works without changing the data."""
        sc = self._sc
        if sc.mask is None:
            sc = sc.with_mask(np.isfinite(sc.unmasked_data[:]))
        return sc

    def _metadata_with_beam(self, beam) -> Metadata:
        """This object's metadata with the beam updated to *beam* (the new resolution)."""
        from dataclasses import replace

        return replace(self.metadata, beam=beam)

    # -- ContextCarrier seam ------------------------------------------------------------
    @property
    def _data_quantity(self) -> u.Quantity:
        return self._sc.filled_data[:]

    def _with_data(self, new_quantity: u.Quantity) -> "Cube":
        sc = _build_cube(new_quantity, self._sc.wcs, self.metadata.beam)
        return Cube(sc, self._metadata_with_unit(new_quantity.unit))

    # -- coordinate-array + validity emission (issue #26) -------------------------------
    @property
    def coordinates(self) -> _coords.AxisCoordinates:
        """Per-axis physical coordinate arrays as ``Quantity``, read from the WCS / spectral
        axis already parsed (no reparse). The absolute frequency is authoritative; the
        per-line Δv is derived from it under the object's convention + rest frequency, and
        raises if that context is absent (no silent guess — ADR-0003)."""
        longitude, latitude = _coords.sky_coordinate_maps(
            self._sc.spatial_coordinate_map
        )
        return _coords.AxisCoordinates(
            frequency=_coords.absolute_frequency(
                self.spectral_axis,
                rest_frequency=self.rest_frequency,
                convention=self.velocity_convention,
            ),
            delta_v=_coords.delta_v(
                self.spectral_axis,
                rest_frequency=self.rest_frequency,
                convention=self.velocity_convention,
            ),
            longitude=longitude,
            latitude=latitude,
            pixel_scale=_coords.pixel_scale(self._sc.wcs),
        )

    @property
    def validity(self) -> _coords.Validity:
        """Validity descriptor: blanked / edge / outside-coverage voxels as ``NaN`` + a
        boolean finite-data mask. Derived from the data, so it travels through a subcube
        slice (mask-of-a-slice == slice-of-the-mask)."""
        return _coords.validity_of(self._data_quantity)

    # -- frequency-authoritative spectral axis + frame transform (issue #24) -------------
    def delta_v_for(self, rest_frequency: u.Quantity) -> u.Quantity:
        """The derived per-channel velocity offset relative to a *chosen* line.

        :attr:`coordinates.delta_v <astrolyze.core._coords.AxisCoordinates.delta_v>` is Δv
        against the primary :attr:`rest_frequency`; broadband work needs Δv per line, so this
        gives Δv against any line's ``rest_frequency`` (e.g. a second species in the band). It
        is derived from the authoritative absolute frequency under the object's convention —
        not assumed — and raises if that convention is absent (ADR-0003)."""
        return _coords.delta_v_relative_to(
            self.coordinates.frequency,
            line_rest_frequency=rest_frequency,
            convention=self.velocity_convention,
        )

    def to_spectral_frame(
        self,
        frame: str,
        *,
        location=None,
        obstime=None,
        target=None,
    ) -> "Cube":
        """Transform the spectral axis to a different reference *frame* (e.g. LSRK->BARYCENT).

        Delegates the maths to astropy :class:`~astropy.coordinates.SpectralCoord` (astrolyze
        stays thin). It is **context-or-raise** (ADR-0003): the observer ``location`` +
        ``obstime`` and the ``target`` sky position are all required — without them a frame
        shift cannot be computed and we raise rather than silently leaving the frame unchanged.
        The current frame (the dataset's ``SPECSYS``) must be known too: an absent
        ``spectral_frame`` is flagged here (lazy enforcement) and raises naming it.

        Returns a new :class:`Cube` on a frequency spectral axis in the requested *frame*, its
        :attr:`~astrolyze.io.Metadata.spectral_frame` updated to *frame* (ADR-0004)."""
        from astropy.wcs import WCS

        shifted_frequency = _coords.to_spectral_frame(
            self.coordinates.frequency,
            frame,
            source_frame=self.metadata.spectral_frame,
            location=location,
            obstime=obstime,
            target=target,
        )
        # Rebuild the spectral WCS on the shifted absolute frequencies; the celestial axes are
        # untouched. spectral-cube's with_spectral_unit keeps freq<->velocity routing through
        # astropy; the frame shift is a (Doppler) rescaling of the linear frequency grid, so we
        # set the reference channel + channel spacing from the shifted axis (CRPIX 1 = chan 0).
        new_sc = self._sc.with_spectral_unit(u.Hz)
        new_wcs = WCS(new_sc.header)
        spec = new_wcs.wcs.spec
        new_wcs.wcs.crpix[spec] = 1.0
        new_wcs.wcs.crval[spec] = shifted_frequency[0].to_value(u.Hz)
        if shifted_frequency.size > 1:
            spacing = (shifted_frequency[1] - shifted_frequency[0]).to_value(u.Hz)
            new_wcs.wcs.cdelt[spec] = spacing
        rebuilt = SpectralCube(data=new_sc.unmasked_data[:], wcs=new_wcs)
        from dataclasses import replace

        _emit("to_spectral_frame", params={"frame": frame})
        return Cube(rebuilt, replace(self.metadata, spectral_frame=frame))

    # -- thin passthroughs --------------------------------------------------------------
    @property
    def spectral_axis(self):
        return self._sc.spectral_axis

    @property
    def unit(self) -> u.UnitBase:
        return self._sc.unit

    @property
    def shape(self) -> tuple[int, ...]:
        return self._sc.shape

    def __repr__(self) -> str:
        obj = self.metadata.object or "?"
        return f"<Cube {obj} {self.shape} [{self.unit}]>"

    # -- spatial tiling / postage-stamp cutouts (issue #30) -----------------------------
    # A cutout iterator: pure SPATIAL slicing that keeps the FULL spectral axis on every
    # tile. astrolyze stays thin (ADR-0004) — each window is just ``self[:, y0:y1, x0:x1]``,
    # so spectral-cube does the slice (and the sub-image WCS) and ``__getitem__`` carries the
    # beam / rest frequency / convention onto the resulting Cube. On a dask-backed cube the
    # slice is a dask graph, so tiling is lazy with no full materialisation (#23). No new
    # resampling/interpolation/padding maths is added here.
    def tiles(self, size, *, stride=None, overlap=0, partial: str = "drop"):
        """Iterate spatial cutout sub-:class:`Cube`\\ s, each keeping the **full** spectral axis.

        A generator over a strided spatial grid — the general postage-stamp / mosaic-tiling op.
        Each yielded cutout is ``self[:, y0:y1, x0:x1]``: a pure spatial slice (no resampling,
        interpolation, or padding) that preserves **all** channels and carries this cube's
        context plus a correct sub-image WCS (ADR-0004). It is **lazy** on a dask-backed cube
        (#23): the cutout is a dask graph, so nothing materialises until you compute it.

        Parameters
        ----------
        size : int or (int, int)
            Spatial tile size in pixels — a square ``n`` or ``(ny, nx)``.
        stride : int or (int, int), optional
            Step between successive tile origins. When given it is authoritative; otherwise it
            is derived from *overlap* as ``size - overlap`` (the natural overlapping-tile
            controls). Must be positive.
        overlap : int or (int, int), default 0
            Convenience for ``stride = size - overlap`` (ignored when *stride* is given);
            ``overlap=0`` is a non-overlapping grid. Must be smaller than *size*.
        partial : {"drop", "keep"}, default "drop"
            The **explicit** edge rule (no silent padding, ADR-0003). ``"drop"`` skips any
            edge window narrower/shorter than *size*; ``"keep"`` yields it as a smaller,
            truncated cutout (still a pure slice — never padded out to *size*).

        Yields
        ------
        Cube
            One spatial cutout per grid position, in row-major (``y`` then ``x``) order.
        """
        ny, nx = _as_pair(size, "size")
        if stride is not None:
            sy, sx = _as_pair(stride, "stride")
        else:
            oy, ox = _as_pair(overlap, "overlap", allow_zero=True)
            sy, sx = ny - oy, nx - ox
        if sy <= 0 or sx <= 0:
            raise ValueError(
                f"tile stride must be positive (got {(sy, sx)}); with overlap, overlap must be "
                "smaller than size so successive tiles advance"
            )
        if partial not in ("drop", "keep"):
            raise ValueError(
                f"partial must be 'drop' or 'keep' (the explicit edge rule), got {partial!r}; "
                "astrolyze never silently pads a partial edge tile (ADR-0003)"
            )

        full_y, full_x = self.shape[1], self.shape[2]
        for y0 in range(0, full_y, sy):
            y1 = min(y0 + ny, full_y)
            if partial == "drop" and (y1 - y0) < ny:
                continue
            for x0 in range(0, full_x, sx):
                x1 = min(x0 + nx, full_x)
                if partial == "drop" and (x1 - x0) < nx:
                    continue
                yield self[:, y0:y1, x0:x1]


def _as_pair(value, name: str, *, allow_zero: bool = False) -> tuple[int, int]:
    """Coerce a tiling parameter to a positive ``(int, int)`` pair (a scalar means square).

    A guard, not silent physics (ADR-0003): a non-integer or non-positive size/stride is a
    programming error and is rejected with a clear message rather than rounded or defaulted."""
    if isinstance(value, (int, np.integer)):
        pair = (int(value), int(value))
    else:
        seq = tuple(value)
        if len(seq) != 2:
            raise ValueError(f"{name} must be an int or a (ny, nx) pair, got {value!r}")
        pair = (int(seq[0]), int(seq[1]))
    floor = 0 if allow_zero else 1
    if pair[0] < floor or pair[1] < floor:
        raise ValueError(
            f"{name} must be {'non-negative' if allow_zero else 'positive'}, got {pair}"
        )
    return pair


def _build_cube(data, wcs, beam, *, unit=None) -> SpectralCube:
    """Construct a SpectralCube, attaching the beam at build time when present.

    A dask-backed array is kept lazy by building a dask cube (``use_dask=True``) and carrying
    its unit via ``meta`` rather than a ``Quantity`` wrap — ``u.Quantity`` would materialise the
    dask array to numpy, defeating the lazy Zarr backend (#23). A unit-carrying ``Quantity``
    (the eager FITS path and the transform methods) is unchanged: ``data`` already states its
    unit, so *unit* is only consulted for the bare-dask case."""
    use_dask = _is_dask(data)
    kwargs = {"wcs": wcs, "use_dask": use_dask}
    if beam is not None:
        kwargs["beam"] = beam
    if use_dask and not hasattr(data, "unit") and unit is not None:
        # Bare dask array: spectral-cube reads the unit from the header meta (a Quantity would
        # have materialised it). dimensionless is recorded as the empty BUNIT it round-trips to.
        bunit = "" if unit == u.dimensionless_unscaled else u.Unit(unit).to_string()
        kwargs["meta"] = {"BUNIT": bunit}
    return SpectralCube(data=data, **kwargs)


def _is_dask(data) -> bool:
    """Whether *data* (possibly a unit-carrying Quantity) wraps a dask array — the lazy path."""
    import dask

    value = getattr(data, "value", data)
    return dask.is_dask_collection(value)


def _is_temperature_unit(unit) -> bool:
    """Whether *unit* is a brightness temperature (K-equivalent) — routes harmonisation."""
    return u.Unit(unit).is_equivalent(u.K)
