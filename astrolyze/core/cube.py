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

from . import _coords
from ._base import ContextCarrier, _emit
from .map import Map
from .spectrum import Spectrum

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
        data = u.Quantity(loaded.data, bunit)
        sc = _build_cube(data, loaded.wcs, loaded.metadata.beam)
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


def _build_cube(data: u.Quantity, wcs, beam) -> SpectralCube:
    """Construct a SpectralCube, attaching the beam at build time when present."""
    if beam is not None:
        return SpectralCube(data=data, wcs=wcs, beam=beam)
    return SpectralCube(data=data, wcs=wcs)
