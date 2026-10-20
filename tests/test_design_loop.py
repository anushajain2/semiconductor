import math

import numpy as np
import pandas as pd
import pytest

from tcid.design_space import DesignSpace
from tcid.search import design_objectives
from tcid.studies import _variance_reduction
from tcid.surrogate import GPSurrogate, aggregate, cross_validate


@pytest.fixture
def space():
    return DesignSpace(
        bounds={"d50": (4e-6, 24e-6), "sigma_ln": (0.05, 0.5), "aspect": (1.0, 4.0), "order": (0.0, 0.95)},
        loading=0.3,
        preferred_angle=math.pi / 2,
    )


def test_unit_round_trip(space):
    for u in space.sobol(16, 1):
        d = space.from_unit(u)
        assert np.allclose(space.to_unit(d), u)
        assert space.contains(d)
        assert d.loading == 0.3


def _synthetic_rows(space, n=24, seeds=3):
    rng = np.random.default_rng(0)
    rows = []
    scen = [("nominal", "nominal", 1e-7, 30.0), ("a", "stress", 1e-8, 20.0), ("b", "stress", 1e-6, 35.0)]
    for i, u in enumerate(space.sobol(n, 3)):
        d = space.from_unit(u)
        for s in range(seeds):
            for name, role, r, kf in scen:
                true = 2.5e-4 * (1 - 0.2 * u[2] * u[3]) * (1 + 0.05 * math.log10(r / 1e-8)) * (30 / kf) ** 0.05
                rows.append({"design_id": f"d{i}", "scenario": name, "scenario_role": role, "seed": s, "status": "ok",
                             "r_layer": true * (1 + 0.01 * rng.standard_normal()), "d50": d.d50, "sigma_ln": d.sigma_ln,
                             "aspect": d.aspect, "order": d.order, "loading_target": 0.3, "r_interface": r,
                             "k_filler_par": kf, "tag": "t", "loading_achieved": 0.3, "geo_order_realized": d.order,
                             "solve_time_s": 1.0})
    return pd.DataFrame(rows)


def test_objectives_and_surrogate(space):
    rows = _synthetic_rows(space)
    obj = design_objectives(rows, n_stress=2)
    assert len(obj) == 24
    assert np.all(obj["robust"] >= obj["mean_stress"])
    agg = aggregate(rows)
    cv = cross_validate(agg, space, lambda: GPSurrogate(n_restarts=0), n_splits=4)
    assert cv["median_rel_error"] < 0.05


def test_variance_reduction():
    df = pd.DataFrame({"a": np.repeat([0, 1], 4), "b": np.tile([0, 1], 4)})
    df["y"] = 3 * df["a"]
    r = _variance_reduction(df, "y", {"a": "a", "b": "b"})
    assert r["variance_reduction_a"] == pytest.approx(1.0)
    assert r["variance_reduction_b"] == pytest.approx(0.0)


def test_gp_noise_is_scaled_to_normalised_units():
    """With normalize_y=True the supplied per-point noise must be divided by var(y)."""
    rng = np.random.default_rng(3)
    X = rng.uniform(size=(30, 2))
    y = 0.2 * np.sin(4 * X[:, 0]) + 0.05 * X[:, 1]
    nv = np.full(30, 1e-4)
    m = GPSurrogate(n_restarts=0).fit(X, y, nv)
    assert np.allclose(m.model.alpha, np.maximum(nv / np.var(y), 1e-8))
    m2 = GPSurrogate(n_restarts=0, normalize_y=False).fit(X, y, nv)
    assert np.allclose(m2.model.alpha, nv)
