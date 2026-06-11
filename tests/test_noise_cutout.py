"""Tests for ``NoiseModel.cutout`` — the spatial cutout of a noise model.

Written first (red/green TDD). Contract: ``model.cutout(y0, y1, x0, x1)`` mirrors
``Cube[:, y0:y1, x0:x1]`` — it slices the spatial σ pattern (and the parent cube) the same way,
keeping the full spectral axis. The motivation is performance: ``.sigma_cube`` reconstructs the σ
field over the WHOLE cube, so reading a small spatial tube's per-voxel σ would otherwise force a
full-cube reconstruction. ``cutout`` slices the model *first* (cheap — only σ_xy / σ_field are
sliced), so ``.sigma_cube`` on the cutout reconstructs just the tube.

Correctness obligation (the numbers must not move): a cutout's σ products equal the full model's
products sliced to the same bbox — same numbers, just cheaper. The spectral σ_v, the scalar, the
ACF and the quality flag are per-cube and unchanged.

Reuses the synthetic fixtures from ``test_noise`` (a stationary/separable noise cube and a
structured/non-separable one) so both representation branches are exercised.
"""

import numpy as np
import astropy.units as u

from astrolyze.core import Cube, Map, NoiseModel

# Reuse the proven fixtures (header schema, separable + non-separable synthetic cubes).
from test_noise import noise_cube, structured_noise_cube  # noqa: F401

# A representative interior bounding box for the 16x16 spatial fixtures (not the full extent).
Y0, Y1, X0, X1 = 3, 11, 5, 13


def _sigma_values(noise_model) -> np.ndarray:
    """The bare σ-cube values (K) of *noise_model* — what the consistency checks compare."""
    return np.asarray(noise_model.sigma_cube._data_quantity.to_value(u.K))


# --------------------------------------------------------------------------------------
# AC: separable consistency — the cutout's σ-cube equals the full σ-cube sliced to the bbox
# --------------------------------------------------------------------------------------
def test_separable_cutout_sigma_cube_matches_full_sliced(noise_cube):
    model = noise_cube.estimate_noise()
    assert model.is_separable is True

    cut = model.cutout(Y0, Y1, X0, X1)
    assert isinstance(cut, NoiseModel)

    full = _sigma_values(model)
    expected = full[:, Y0:Y1, X0:X1]
    got = _sigma_values(cut)

    # Same shape as the bbox (full spectral axis kept) and identical numbers — just cheaper.
    assert got.shape == (model.sigma_cube.shape[0], Y1 - Y0, X1 - X0)
    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=1e-12)


# --------------------------------------------------------------------------------------
# AC: per-cube quantities are unchanged; only the spatial products take the bbox shape
# --------------------------------------------------------------------------------------
def test_separable_cutout_keeps_per_cube_quantities(noise_cube):
    model = noise_cube.estimate_noise()
    cut = model.cutout(Y0, Y1, X0, X1)

    # The spatial σ map takes the bbox spatial shape ...
    assert isinstance(cut.sigma_map, Map)
    assert cut.sigma_map.data.shape == (Y1 - Y0, X1 - X0)

    # ... while σ_v, the ACF, the scalar and the quality are per-cube and unchanged.
    np.testing.assert_allclose(
        np.asarray(cut.sigma_spectrum.flux.to_value(u.K)),
        np.asarray(model.sigma_spectrum.flux.to_value(u.K)),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.asarray(cut.spectral_acf.flux.value),
        np.asarray(model.spectral_acf.flux.value),
        rtol=1e-12,
        atol=1e-12,
    )
    assert cut.scalar.to_value(u.K) == model.scalar.to_value(u.K)
    assert cut.quality is model.quality
    # Provenance (method / version) rides along too.
    assert cut.method == model.method
    assert cut.version == model.version


# --------------------------------------------------------------------------------------
# AC: full (non-separable) consistency — same equality, on the dense σ_field branch
# --------------------------------------------------------------------------------------
def test_full_cutout_sigma_cube_matches_full_sliced(structured_noise_cube):
    model = structured_noise_cube.estimate_noise()
    assert model.is_separable is False  # this fixture trips the full-3D fallback

    cut = model.cutout(Y0, Y1, X0, X1)
    assert cut.is_separable is False

    full = _sigma_values(model)
    expected = full[:, Y0:Y1, X0:X1]
    got = _sigma_values(cut)

    assert got.shape == (model.sigma_cube.shape[0], Y1 - Y0, X1 - X0)
    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=1e-12)


# --------------------------------------------------------------------------------------
# AC: identity — a full-extent cutout reproduces the original σ-cube exactly
# --------------------------------------------------------------------------------------
def test_full_extent_cutout_is_identity(noise_cube, structured_noise_cube):
    for cube in (noise_cube, structured_noise_cube):
        model = cube.estimate_noise()
        _, ny, nx = model.sigma_cube.shape
        cut = model.cutout(0, ny, 0, nx)
        np.testing.assert_allclose(
            _sigma_values(cut), _sigma_values(model), rtol=1e-12, atol=1e-12
        )


# --------------------------------------------------------------------------------------
# AC: the cutout's parent cube is the spatially-sliced cube (sub-WCS/shape carried)
# --------------------------------------------------------------------------------------
def test_cutout_parent_cube_is_the_sliced_cube(noise_cube):
    model = noise_cube.estimate_noise()
    cut = model.cutout(Y0, Y1, X0, X1)
    assert isinstance(cut._cube, Cube)
    assert cut._cube.shape == (model._cube.shape[0], Y1 - Y0, X1 - X0)
