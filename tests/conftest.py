import pytest

from tcid.config import Design, GeometrySettings, Layer, MeshSettings
from tcid.geometry import generate
from tcid.mesh import mesh_particles


@pytest.fixture(scope="session")
def small_layer():
    return Layer(thickness=50e-6, width=100e-6, margin=3e-6)


@pytest.fixture(scope="session")
def small_geometry(small_layer):
    d = Design(loading=0.30, d50=8e-6, sigma_ln=0.25, aspect=2.0, order=0.5)
    g = generate(d, small_layer, GeometrySettings(), seed=5)
    assert g.ok, g.reason
    return g


@pytest.fixture(scope="session")
def small_mesh(small_geometry):
    return mesh_particles(small_geometry, MeshSettings(size_factor=2.0))
