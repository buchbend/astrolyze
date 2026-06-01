#!/usr/bin/env python
"""The astrolyze tracer bullet, as a plain script (PRD 0001 acceptance).

This is the whole vertical slice — load a PPV cube, collapse it to a velocity-integrated
map, put that map in the units you want, and display it the house way — written exactly as
a human would write it. There is no AI-only path: a coding agent drives astrolyze through
this same ordinary, reviewable Python (ADR-0011/0012). The corresponding one-liner from the
shell is ``astrolyze moment0 <cube.fits>``.

Run it on the full NGC 628 CO(2-1) cube::

    python examples/tracer_ngc628.py /path/to/ngc0628_12m+7m+tp_co21.fits
    # or point $ASTROLYZE_TRACER_CUBE at the cube and run with no argument

It writes ``ngc0628_co21_moment0.png`` (override with a second argument).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from astrolyze.core import Cube
from astrolyze.io import load


def main(cube_path: str, out_path: str) -> None:
    # 1. Load: the FITS header is authoritative. astrolyze reads the schema (beam, rest
    #    frequency, and — from CTYPE3=VRAD — the radio velocity convention) and the load is
    #    lazy, so the heavy array is only touched when we compute on it.
    cube = Cube.from_loaded(load(cube_path))
    print(f"loaded {cube!r}")  # e.g. <Cube NGC0628 (98, 1600, 1600) [K]>

    # 2. Velocity-integrate to moment 0 -> a Map. The beam / rest frequency / convention ride
    #    along on the result; we never re-state them. For a brightness-temperature (K) cube
    #    this gives K·(m/s).
    integrated = cube.moment0()

    # 3. Put it in the units we want. The Map supplies its own context to the unit hub, so
    #    K·(m/s) -> K·(km/s) is a clean rescale with nothing to spell out. (A genuinely
    #    ambiguous brightness conversion would require an explicit Rayleigh-Jeans/Planck
    #    scale here — astrolyze refuses to guess it; ADR-0003.)
    integrated = integrated.to("K km/s")
    print(f"moment-0 map: {integrated.unit}")

    # 4. Display the house way: cividis on WCS axes, the beam drawn from the object, a
    #    colorbar labelled in the data's units. Plotting never touches your global matplotlib
    #    settings (the style is applied locally; ADR-0005).
    fig, _ = integrated.plot()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ASTROLYZE_TRACER_CUBE")
    if not path:
        raise SystemExit(
            "usage: python examples/tracer_ngc628.py <cube.fits> [out.png]\n"
            "   (or set $ASTROLYZE_TRACER_CUBE to the cube path)"
        )
    out = sys.argv[2] if len(sys.argv) > 2 else f"{Path(path).stem}_moment0.png"
    main(path, out)
