"""Smoke tests for ``astrolyze collection list PATH`` (issue #57).

The CLI mirrors the Python API: ``astrolyze collection list`` opens a published corpus through
the public :class:`~astrolyze.collection.Collection` and prints the object-first overview as a
rich table — the shell path to the same list a human reaches in Python. As elsewhere these are
contract smoke tests through the typer runner (no browser, no large files); the synthetic corpus
fixture is shared with the library suite.
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


@pytest.fixture
def corpus(tmp_path):
    root = tmp_path / "ism_corpus"
    root.mkdir()
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"], h["TELESCOP"], h["BUNIT"] = "NGC3521", "IRAM30M", "K"
    h["RESTFRQ"] = (REST_CO21, "Hz")
    beam = radio_beam.Beam(major=13 * u.arcsec, minor=11 * u.arcsec, pa=30 * u.deg)
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO"
    data = np.arange(3 * 4 * 5, dtype="float32").reshape(3, 4, 5)
    (root / "l1").mkdir()
    Cube.from_loaded(_loaded(data, h)).to_zarr(root / "l1")
    store = next((root / "l1").glob("*.zarr"))

    from astrolyze.collection.catalog import CATALOG_COLUMNS

    row = {
        "object": "NGC3521",
        "survey": "HERACLES",
        "telescope": "IRAM30M",
        "species": "CO",
        "transition": "2-1",
        "rest_frequency_hz": REST_CO21,
        "beam_major_arcsec": 13.0,
        "beam_minor_arcsec": 11.0,
        "beam_pa_deg": 30.0,
        "bunit": "K",
        "store_path": store.relative_to(root).as_posix(),
        "content_checksum": "sha256:x",
        "ra_deg": 170.0,
        "dec_deg": -0.04,
        "radius_deg": 0.01,
        "catalog_schema_version": "1.0",
    }
    table = pa.table({k: [row.get(k)] for k in CATALOG_COLUMNS})
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
