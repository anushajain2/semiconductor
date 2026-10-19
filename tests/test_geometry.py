import math

import numpy as np
import pytest

from tcid.config import Design, GeometrySettings, Layer
from tcid.geometry import generate, kappa_from_order, realized_order, sample_orientations, validate


def test_order_parameter_inversion():
    rng = np.random.default_rng(0)
    for s in (0.2, 0.6, 0.9):
        k = kappa_from_order(s)
        th = sample_orientations(rng, 20000, s, math.pi / 2)
        assert realized_order(th, math.pi / 2) == pytest.approx(s, abs=0.02)
        assert k > 0


def test_random_orientation_has_zero_order():
    th = sample_orientations(np.random.default_rng(1), 20000, 0.0, math.pi / 2)
    assert abs(realized_order(th, math.pi / 2)) < 0.03


def test_generated_geometry_is_valid_and_recorded(small_geometry):
    v = validate(small_geometry)
    assert v["valid"]
    assert v["min_gap_lower_bound"] >= small_geometry.settings["min_gap"] * 0.995
    assert abs(small_geometry.stats["loading_ellipse"] - 0.30) <= 0.005
    assert small_geometry.stats["n_particles"] == small_geometry.n


def test_reproducible_with_seed(small_layer):
    d = Design(loading=0.25, d50=6e-6, sigma_ln=0.3, aspect=1.5, order=0.3)
    g1 = generate(d, small_layer, GeometrySettings(), 3)
    g2 = generate(d, small_layer, GeometrySettings(), 3)
    g3 = generate(d, small_layer, GeometrySettings(), 4)
    assert np.array_equal(g1.cx, g2.cx)
    assert not np.array_equal(g1.cx, g3.cx)


def test_impossible_requests_are_rejected_not_altered(small_layer):
    # particles longer than the layer
    g = generate(Design(0.3, 40e-6, 0.1, 4.0, 0.9), small_layer, GeometrySettings(), 0)
    assert not g.ok and "truncated" in g.reason
    # unpackable loading
    g = generate(Design(0.75, 6e-6, 0.1, 1.0, 0.0), small_layer, GeometrySettings(max_attempts=300, restarts=1), 0)
    assert not g.ok and "packing infeasible" in g.reason


def test_particles_respect_margins(small_geometry):
    lay = small_geometry.layer
    from tcid.geometry import half_extents

    ex, ey = half_extents(small_geometry.a, small_geometry.b, small_geometry.theta)
    assert np.all(small_geometry.cy - ey >= lay["margin"] - 1e-12)
    assert np.all(small_geometry.cy + ey <= lay["thickness"] - lay["margin"] + 1e-12)
