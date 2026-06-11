"""Smoke tests for ``astrolyze collection list|describe|query`` (issues #57, #60).

The CLI mirrors the Python API: ``astrolyze collection list`` opens a published corpus through
the public :class:`~astrolyze.collection.Collection` and prints the object-first overview as a
rich table — the shell path to the same list a human reaches in Python. ``describe OBJECT PATH``
and ``query PATH --species … --survey …`` mirror the #60 drill-down + filter operations. As
elsewhere these are contract smoke tests through the typer runner (no browser, no large files);
the synthetic corpus fixture is shared with the library suite.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits
from typer.testing import CliRunner

from astrolyze.cli import app
from astrolyze.core import Cube

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

runner = CliRunner()

REST_CO21 = 230.538e9


def _all_output(result) -> str:
    out = result.stdout or ""
    try:
        out += result.stderr or ""
    except ValueError:
        pass
    return out


def _loaded(data, header):
    from astropy.wcs import WCS

    from astrolyze.io.access import LoadedData
    from astrolyze.io.schema import Metadata

    return LoadedData(
        data=data,
        wcs=WCS(header),
        metadata=Metadata.from_header(header),
        path="<memory>",
        header_string=header.tostring(),
        header=header,
    )


def _header(*, obj, telescope, species, rest_hz, bmaj):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"], h["TELESCOP"], h["BUNIT"] = obj, telescope, "K"
    h["RESTFRQ"] = (rest_hz, "Hz")
    beam = radio_beam.Beam(
        major=bmaj * u.arcsec, minor=(bmaj - 2) * u.arcsec, pa=30 * u.deg
    )
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = species
    return h, beam


@pytest.fixture
def corpus(tmp_path):
    """Two stores: NGC3521 in HERACLES (CO 2-1) and NGC0628 in THINGS (HI 21cm)."""
    root = tmp_path / "ism_corpus"
    root.mkdir()
    (root / "l1").mkdir()
    REST_HI = 1.420405751e9

    specs = [
        dict(
            obj="NGC3521",
            survey="HERACLES",
            telescope="IRAM30M",
            species="CO",
            transition="2-1",
            rest_hz=REST_CO21,
            bmaj=13.0,
        ),
        dict(
            obj="NGC0628",
            survey="THINGS",
            telescope="VLA",
            species="HI",
            transition="21cm",
            rest_hz=REST_HI,
            bmaj=11.0,
        ),
    ]

    from astrolyze.collection.catalog import CATALOG_COLUMNS

    rows = []
    for spec in specs:
        h, beam = _header(
            obj=spec["obj"],
            telescope=spec["telescope"],
            species=spec["species"],
            rest_hz=spec["rest_hz"],
            bmaj=spec["bmaj"],
        )
        data = np.arange(3 * 4 * 5, dtype="float32").reshape(3, 4, 5)
        sub = root / "l1" / spec["survey"].lower()
        sub.mkdir()
        Cube.from_loaded(_loaded(data, h)).to_zarr(sub)
        store = next(sub.glob("*.zarr"))
        rows.append(
            {
                "object": spec["obj"],
                "survey": spec["survey"],
                "telescope": spec["telescope"],
                "species": spec["species"],
                "transition": spec["transition"],
                "rest_frequency_hz": spec["rest_hz"],
                "beam_major_arcsec": beam.major.to_value(u.arcsec),
                "beam_minor_arcsec": beam.minor.to_value(u.arcsec),
                "beam_pa_deg": beam.pa.to_value(u.deg),
                "bunit": "K",
                "store_path": store.relative_to(root).as_posix(),
                "content_checksum": "sha256:x",
                "ra_deg": 170.0,
                "dec_deg": -0.04,
                "radius_deg": 0.01,
                "catalog_schema_version": "1.0",
            }
        )

    table = pa.table({k: [row.get(k) for row in rows] for k in CATALOG_COLUMNS})
    table = table.replace_schema_metadata({"catalog_schema_version": "1.0"})
    pq.write_table(table, root / "catalog.parquet")
    return root


def test_collection_list_prints_object_first_table(corpus):
    result = runner.invoke(app, ["collection", "list", str(corpus)])
    assert result.exit_code == 0, _all_output(result)
    out = _all_output(result)
    assert "NGC3521" in out
    assert "HERACLES" in out
    assert "CO" in out


def test_collection_list_catalog_less_empty_directory_reports_no_datasets(tmp_path):
    """Since #61, a catalog-less directory is *scanned* rather than rejected: an empty one has no
    stores to find, so the CLI reports "no datasets" and exits cleanly (the scan fallback, not an
    error). A directory with stores but no catalog.parquet lists the scanned cubes instead."""
    empty = tmp_path / "empty"
    empty.mkdir()
    result = runner.invoke(app, ["collection", "list", str(empty)])
    assert result.exit_code == 0, _all_output(result)
    assert "no datasets" in _all_output(result).lower()


# --------------------------------------------------------------------------------------
# collection describe OBJECT PATH (issue #60)
# --------------------------------------------------------------------------------------
def test_collection_describe_renders_detail_table(corpus):
    result = runner.invoke(app, ["collection", "describe", "NGC3521", str(corpus)])
    assert result.exit_code == 0, _all_output(result)
    out = _all_output(result)
    assert "NGC3521" in out
    assert "HERACLES" in out
    assert "CO" in out
    assert "IRAM30M" in out


def test_collection_describe_unknown_object_exits_nonzero(corpus):
    result = runner.invoke(app, ["collection", "describe", "NGC9999", str(corpus)])
    assert result.exit_code == 1, _all_output(result)
    assert "NGC9999" in _all_output(result)


def test_collection_describe_deep_shows_velocity_columns(corpus):
    result = runner.invoke(
        app, ["collection", "describe", "NGC3521", str(corpus), "--deep"]
    )
    assert result.exit_code == 0, _all_output(result)
    # The deep-only channel width (2 km/s on the synthetic axis) appears.
    assert "2.0" in _all_output(result) or "2.00" in _all_output(result)


# --------------------------------------------------------------------------------------
# collection query PATH --species/--survey/... (issue #60)
# --------------------------------------------------------------------------------------
def test_collection_query_filters_by_species(corpus):
    result = runner.invoke(app, ["collection", "query", str(corpus), "--species", "HI"])
    assert result.exit_code == 0, _all_output(result)
    out = _all_output(result)
    assert "NGC0628" in out
    assert "THINGS" in out
    # The CO/HERACLES store is filtered out.
    assert "HERACLES" not in out


def test_collection_query_no_match_reports_no_records(corpus):
    result = runner.invoke(
        app, ["collection", "query", str(corpus), "--species", "SiO"]
    )
    assert result.exit_code == 0, _all_output(result)
    assert "no" in _all_output(result).lower()


def test_collection_query_multi_filter(corpus):
    result = runner.invoke(
        app,
        ["collection", "query", str(corpus), "--species", "CO", "--survey", "HERACLES"],
    )
    assert result.exit_code == 0, _all_output(result)
    assert "NGC3521" in _all_output(result)
