"""Regression tests: the recommendation interface must compare physical values, not names."""

from pathlib import Path

import pytest
import yaml

from tcid.recommend import recommend

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "phase1"
pytestmark = pytest.mark.skipif(not (RESULTS / "final" / "final_rows.parquet").exists(),
                                reason="stored study results not present")


@pytest.fixture
def scen():
    with open(ROOT / "configs" / "scenarios_phase1.yaml") as fh:
        return yaml.safe_load(fh)


def _ask(tmp_path, scenarios: dict, thickness="100 um", loading=0.30):
    p = tmp_path / "scen.yaml"
    p.write_text(yaml.safe_dump(scenarios))
    return recommend(str(RESULTS), "alumina_epoxy_dev_2d", thickness, loading, str(p), str(tmp_path / "out"),
                     config=str(ROOT / "configs" / "study_phase1.yaml"))


def test_unchanged_request_is_supported(tmp_path, scen, monkeypatch):
    monkeypatch.chdir(ROOT)
    assert _ask(tmp_path, scen)["supported"]


def test_changed_interface_resistance_with_same_names_is_rejected(tmp_path, scen, monkeypatch):
    # reviewer's case: every interface resistance moved to 1e-3 m2K/W, scenario names kept
    monkeypatch.chdir(ROOT)
    for s in scen["scenarios"]:
        s["r_interface"] = "1.0e-3 m2K/W"
    res = _ask(tmp_path, scen)
    assert not res["supported"]
    assert any("r_interface" in p for p in res["problems"])


def test_changed_filler_conductivity_is_rejected(tmp_path, scen, monkeypatch):
    monkeypatch.chdir(ROOT)
    scen["scenarios"][0]["k_filler"] = "31 W/mK"
    assert not _ask(tmp_path, scen)["supported"]


def test_renamed_scenario_with_identical_values_is_supported(tmp_path, scen, monkeypatch):
    monkeypatch.chdir(ROOT)
    scen["scenarios"][0]["name"] = "renamed"
    assert _ask(tmp_path, scen)["supported"]


def test_changed_thickness_is_rejected(tmp_path, scen, monkeypatch):
    monkeypatch.chdir(ROOT)
    assert not _ask(tmp_path, scen, thickness="120 um")["supported"]
