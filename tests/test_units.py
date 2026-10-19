import math

import pytest

from tcid.units import from_si, to_si, volume_to_weight_fraction, weight_to_volume_fraction


def test_lengths_and_resistances():
    assert to_si("100 um") == pytest.approx(100e-6)
    assert to_si("0.05 cm2K/W") == pytest.approx(5e-6)
    assert to_si("3 mm2K/W") == pytest.approx(3e-6)
    assert to_si(2.0, "mm") == pytest.approx(2e-3)
    assert to_si({"value": 90, "unit": "deg"}) == pytest.approx(math.pi / 2)
    assert from_si(1e-6, "um") == pytest.approx(1.0)


def test_unknown_unit():
    with pytest.raises(ValueError):
        to_si("3 furlongs")


def test_fraction_round_trip():
    phi = weight_to_volume_fraction(0.5, 2200.0, 1200.0)
    assert 0 < phi < 0.5  # denser filler -> lower volume fraction than weight fraction
    assert volume_to_weight_fraction(phi, 2200.0, 1200.0) == pytest.approx(0.5)
