# Test data provenance

## `ngc0628_co21_cutout.fits.gz`

A small **128 × 128 pixel × 50 channel** cutout of the PHANGS–ALMA CO(2–1) data cube of
the nearby spiral galaxy **NGC 628**, taken over the inner-galaxy region that carries the
most integrated line flux. It is committed solely as a **test fixture** — it lets the
tracer-path smoke test (`tests/test_tracer.py::test_tracer_real_cutout`) run end-to-end on
*real* data in CI without shipping the full ~958 MB cube.

- **Source:** the public PHANGS–ALMA combined `12m+7m+tp` CO(2–1) cube for NGC 628
  (`ngc0628_12m+7m+tp_co21.fits`), part of the PHANGS–ALMA survey.
- **How it was made:** `python tests/data/make_cutout.py <full_cube.fits>` — see that script
  for the exact spatial/spectral window. The slice preserves the full header schema (beam
  `BMAJ/BMIN/BPA`, `RESTFRQ`, `CTYPE3=VRAD`, `SPECSYS`, `BUNIT=K`) and shifts the WCS
  reference pixels accordingly, so astrolyze parses it as *complete*.
- **No reprocessing:** voxel values are the original calibrated brightness temperatures
  (K); only a spatial+spectral crop and gzip were applied.

### Attribution

PHANGS–ALMA data are public. If you use the science content, please cite the survey —
Leroy et al. 2021, ApJS, 257, 43 (PHANGS–ALMA: survey description) — and include the
standard ALMA acknowledgement. This cutout is a derived, downsampled product used here for
software testing only and is not intended for scientific use.
