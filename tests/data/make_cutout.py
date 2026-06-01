"""Generate the committed NGC 628 test cutout from the full PHANGS-ALMA cube.

The full cube (1600x1600x98, ~958 MB) is far too large for a public repo, so the
real-data smoke test ships a small spatial+spectral cutout instead. This script
documents *exactly* how that fixture was produced so it is reproducible and the
provenance is auditable (see PROVENANCE.md).

Run from the repo root with the full cube on disk::

    python tests/data/make_cutout.py /path/to/ngc0628_12m+7m+tp_co21.fits

It writes ``tests/data/ngc0628_co21_cutout.fits.gz``. The cutout is taken over the
128x128 spatial window with the most integrated CO flux and the channels carrying the
line, and it preserves the full header schema (beam, RESTFRQ, CTYPE3=VRAD, SPECSYS,
BUNIT) so astrolyze parses it as *complete* and the tracer runs without raising.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

# Cutout geometry (fixed so the fixture is deterministic / reproducible).
SPATIAL = 128  # cutout is SPATIAL x SPATIAL pixels
ZLO, ZHI = 30, 80  # channels [ZLO, ZHI) bracketing the CO(2-1) line (peaks at ch 52)
STRIDE = 32  # coarse grid for the max-flux window search

OUT = Path(__file__).resolve().parent / "ngc0628_co21_cutout.fits.gz"


def _best_window(line_flux: np.ndarray, size: int, stride: int) -> tuple[int, int]:
    """Top-left (y0, x0) of the size x size window maximising integrated flux."""
    ny, nx = line_flux.shape
    best, best_yx = -np.inf, (0, 0)
    cumulative = np.nan_to_num(line_flux, nan=0.0)
    for y0 in range(0, ny - size + 1, stride):
        for x0 in range(0, nx - size + 1, stride):
            total = cumulative[y0 : y0 + size, x0 : x0 + size].sum()
            if total > best:
                best, best_yx = total, (y0, x0)
    return best_yx


def main(src: str) -> None:
    with fits.open(src, memmap=True) as hdul:
        hdu = next(h for h in hdul if h.data is not None)
        header = hdu.header.copy()
        # Integrated line emission over the spectral window picks the spatial window.
        line = np.nansum(hdu.data[ZLO:ZHI], axis=0)
        y0, x0 = _best_window(line, SPATIAL, STRIDE)
        cut = np.array(hdu.data[ZLO:ZHI, y0 : y0 + SPATIAL, x0 : x0 + SPATIAL])

    # Slice the WCS by shifting the reference pixels (FITS CRPIX is 1-based, the slice
    # offsets are 0-based, so CRPIX_new = CRPIX_old - offset).
    header["NAXIS1"] = SPATIAL
    header["NAXIS2"] = SPATIAL
    header["NAXIS3"] = ZHI - ZLO
    header["CRPIX1"] = header["CRPIX1"] - x0
    header["CRPIX2"] = header["CRPIX2"] - y0
    header["CRPIX3"] = header["CRPIX3"] - ZLO

    fits.writeto(OUT, cut.astype(np.float32), header, overwrite=True)
    size_mb = OUT.stat().st_size / 1e6
    print(f"wrote {OUT} ({cut.shape}, {size_mb:.2f} MB) from window y0={y0} x0={x0}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python tests/data/make_cutout.py <full_cube.fits>")
    main(sys.argv[1])
