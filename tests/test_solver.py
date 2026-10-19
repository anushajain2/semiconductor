import math

import numpy as np
import pytest

from tcid import analytic
from tcid.mesh import FILLER, MATRIX, mesh_particles, mesh_strips
from tcid.solver import LayerModel, MaterialParams

KM, KF = 0.2, 30.0


def test_homogeneous_layer_exact():
    m = mesh_strips([0, 80e-6], [MATRIX], 160e-6, 10e-6)
    s = LayerModel(m).solve(MaterialParams(KM, KF, KF, 0.0))
    assert s.r_layer == pytest.approx(80e-6 / KM, rel=5e-3)
    assert s.heat_balance_error < 1e-3


@pytest.mark.parametrize("r", [0.0, 1e-7, 1e-6, 1e-5])
def test_layered_series_with_interface_resistance(r):
    lv = [0, 30e-6, 60e-6, 100e-6]
    m = mesh_strips(lv, [MATRIX, FILLER, MATRIX], 200e-6, 5e-6)
    s = LayerModel(m).solve(MaterialParams(KM, KF, KF, r))
    exact = analytic.series_resistance(np.diff(lv), [KM, KF, KM], [r, r])
    assert s.r_layer == pytest.approx(exact, rel=1e-2)


def test_parallel_strips_interfaces_carry_no_flux():
    xl = [0, 40e-6, 100e-6]
    m = mesh_strips(xl, [MATRIX, FILLER], 50e-6, 5e-6, orientation="parallel")
    for r in (0.0, 1e-5):
        s = LayerModel(m).solve(MaterialParams(KM, KF, KF, r))
        assert s.r_layer == pytest.approx(analytic.parallel_resistance(np.diff(xl), [KM, KF], 50e-6), rel=1e-2)


def test_robin_contacts_add_in_series():
    m = mesh_strips([0, 100e-6], [MATRIX], 200e-6, 10e-6)
    s = LayerModel(m).solve(MaterialParams(KM, KF, KF, 0.0, r_hot=2e-5, r_cold=5e-6))
    assert s.r_assembly == pytest.approx(100e-6 / KM + 2.5e-5, rel=1e-6)


def test_perfect_contact_requires_tied_dofs(small_mesh):
    lm = LayerModel(small_mesh)
    split = lm.disc(True)
    tied = lm.disc(False)
    assert split.n_dof > tied.n_dof
    assert split.interface_length == pytest.approx(tied.interface_length)


def test_energy_balance_and_monotonicity(small_mesh):
    lm = LayerModel(small_mesh)
    vals = []
    for r in (0.0, 1e-8, 1e-6, 1e-5, 1e-3):
        s = lm.solve(MaterialParams(KM, KF, KF, r))
        assert s.heat_balance_error < 1e-3
        assert s.flux_consistency_error < 1e-6
        vals.append(s.r_layer)
    assert np.all(np.diff(vals) > 0)


def test_pure_phase_limit(small_mesh):
    s = LayerModel(small_mesh).solve(MaterialParams(KM, KM, KM, 0.0))
    assert s.r_layer == pytest.approx(small_mesh.flow_length / KM, rel=1e-9)


def test_wiener_bounds(small_mesh):
    s = LayerModel(small_mesh).solve(MaterialParams(KM, KF, KF, 0.0))
    k = small_mesh.flow_length / s.r_layer
    lo, hi = analytic.wiener_bounds(small_mesh.filler_fraction(), KM, KF)
    assert lo <= k <= hi


def test_rigid_rotation_with_anisotropic_tensor(small_mesh):
    mp = MaterialParams(KM, 300.0, 5.0, 1e-7)
    r0 = LayerModel(small_mesh).solve(mp).r_layer
    for ang in (0.3, 1.2, math.pi / 2, -2.0):
        assert LayerModel(small_mesh.rotated(ang)).solve(mp).r_layer == pytest.approx(r0, rel=1e-9)
    bad = small_mesh.rotated(math.pi / 2)
    bad.region_theta = small_mesh.region_theta.copy()
    assert abs(LayerModel(bad).solve(mp).r_layer - r0) / r0 > 1e-3


def test_rotated_geometry_remeshed(small_geometry, small_mesh):
    from tcid.config import MeshSettings

    mp = MaterialParams(KM, 300.0, 5.0, 1e-7)
    r0 = LayerModel(small_mesh).solve(mp).r_layer
    m90 = mesh_particles(small_geometry.rotated90(), MeshSettings(size_factor=2.0), flow_axis="x")
    assert LayerModel(m90).solve(mp).r_layer == pytest.approx(r0, rel=1e-2)


def test_no_unresolved_bridges(small_mesh):
    assert small_mesh.info["bridge_elements"] == 0
    assert small_mesh.info["margin_bridges"] == 0
    assert small_mesh.info["mesh_ok"]
