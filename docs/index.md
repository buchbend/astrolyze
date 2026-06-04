# astrolyze documentation

astrolyze is a thin, opinionated layer over the modern Astropy ecosystem
(`astropy`, `spectral-cube`, `radio_beam`, `reproject`, `specutils`) for radio/sub-mm
position–position–velocity (PPV) cubes and spectra. It exists to make data handling, unit
conversion, coordinates, and I/O *consistent and correct* — never silently guessing the
physics that quietly biases this science.

This documentation is **task- and concept-oriented**: each page below tells you what a
capability is, *why* it works the way it does, and shows a short snippet you can run. It is
not the API reference (the docstrings are) and not a first-run walkthrough (that is the
[tutorial notebook](../examples/tutorial.ipynb)).

## Where to start

- **New here?** Read the [README quickstart](../README.md#quickstart), then open the
  interactive [`examples/tutorial.ipynb`](../examples/tutorial.ipynb) (`pip install -e ".[notebook]"`).
- **Want to know what a feature does and why?** Use the capability guides below.
- **Want the decisions behind the design?** See the [Architecture Decision Records](adr/).

## Capability guides

| Guide | What it covers | Embodies |
| --- | --- | --- |
| [I/O and metadata](guide/io-and-metadata.md) | Loading FITS lazily; the typed `Metadata` schema and its two projections (FITS header and a JSON-able `attrs` dict); the backend-neutral `LoadedData` seam | ADR-0006, ADR-0003 |
| [Coordinates and validity](guide/coordinates-and-validity.md) | Per-axis physical coordinate arrays (`Cube.coordinates`) and the validity descriptor (`Cube.validity`), surfaced read-only from the WCS | ADR-0004, ADR-0003 |
| [No silent physics](guide/no-silent-physics.md) | The typed-insufficiency guard: `MissingContextError`, the structured `Insufficiency`/`ContextGap` descriptor, and the non-raising probe `Cube.can_convert_to(...)` | ADR-0003, ADR-0013 |
| [Beam and channel matching](guide/beam-and-channel-matching.md) | Smoothing to a coarser resolution: `convolve_to_beam`, `spectral_bin`, `spectral_smooth_to`, `match_to`, and the lossy-direction guard | ADR-0003, ADR-0004 |

## How the pieces fit

```
load(path)  ->  LoadedData (data + WCS + Metadata + verbatim header string)
                     |
        Cube.from_loaded(...)  /  Map.from_loaded(...)  /  Spectrum.from_oned(...)
                     |
   context-carrying objects: .coordinates  .validity  .to(unit)  .can_convert_to(unit)
                             convolve_to_beam / spectral_bin / spectral_smooth_to / match_to
```

The single source of truth for a dataset's physical context (rest frequency, velocity
convention, beam, brightness unit, …) is the `Metadata` dataclass. The FITS header and the
`attrs` dict are two *projections* of it; the core objects (`Cube`/`Map`/`Spectrum`) carry it
and never re-parse the header.

## Architecture decisions

The [`adr/`](adr/) directory holds the Architecture Decision Records that govern the public
toolkit's design. The guides above cite them by number. The ADRs are mirrored from a private
canonical repository and are read-only here — see [`adr/README.md`](adr/README.md) for the
boundary and why the numbering has gaps.
