import math

import numpy as np
import pytest
from scipy import stats as sps

from tcid.stats import bootstrap, first_crossing, holm, welch


def test_welch_matches_scipy():
    rng = np.random.default_rng(0)
    a, b = rng.normal(1.0, 0.1, 6), rng.normal(1.15, 0.2, 9)
    t = welch(a, b)
    ref = sps.ttest_ind(a, b, equal_var=False)
    assert t["p"] == pytest.approx(ref.pvalue, rel=1e-9)
    assert t["ci_low"] < t["diff"] < t["ci_high"]


def test_holm_is_monotone_and_bounded():
    adj = holm({"a": 0.01, "b": 0.04, "c": 0.03})
    assert adj == pytest.approx({"a": 0.03, "c": 0.06, "b": 0.06})
    assert all(0 <= v <= 1 for v in adj.values())


def test_bootstrap_covers_true_difference():
    rng = np.random.default_rng(1)
    a, b = rng.normal(0.0, 1.0, 40), rng.normal(1.0, 1.0, 40)
    r = bootstrap({"a": a, "b": b}, lambda a, b: b - a, n=2000)
    assert r["ci_low"] < 1.0 < r["ci_high"]


def test_first_crossing_log_interpolation():
    x = np.logspace(-7, -5, 3)  # 1e-7, 1e-6, 1e-5
    assert first_crossing(x, np.array([-1.0, 1.0, 2.0])) == pytest.approx(10 ** -6.5)
    assert math.isnan(first_crossing(x, np.array([1.0, 2.0, 3.0])))
