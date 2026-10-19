"""Command-line interface.

    tcid verify                      solver acceptance gates
    tcid geometry --d50 10um ...     generate (and plot) one realization
    tcid solve case.yaml             solve one configured case
    tcid pilot | initial | surrogate | search | final | questions | gate
    tcid figures | report            reproduce figures and the technical report
    tcid recommend ...               minimal user workflow
    tcid all                         run the complete Phase 1 pipeline
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml


def _study(a):
    from .studies import Study

    return Study(a.config, a.out, a.cache, a.workers)


def cmd_verify(a):
    from .verification import run_all

    res = run_all(Path(a.out).parent / "verification")
    ok = all(r["passed"] for r in res)
    print("all solver gates passed" if ok else "SOME SOLVER GATES FAILED")
    return 0 if ok else 1


def cmd_geometry(a):
    from .config import Design, GeometrySettings, Layer
    from .geometry import generate, validate
    from .units import to_si

    d = Design(loading=a.loading, d50=to_si(a.d50), sigma_ln=a.sigma, aspect=a.aspect, order=a.order)
    g = generate(d, Layer(thickness=to_si(a.thickness), width=to_si(a.width)), GeometrySettings(), a.seed)
    print(json.dumps({"status": g.status, "reason": g.reason, **g.stats, **(validate(g) if g.ok else {})}, indent=2, default=float))
    if g.ok and a.plot:
        from .figures import plot_geometry

        plot_geometry(g, a.plot)
        print(f"wrote {a.plot}")
    return 0 if g.ok else 2


def cmd_solve(a):
    """Solve one case described by a YAML file (design, seed, optional layer/mesh/scenarios)."""
    from .config import design_from_dict, geometry_from_dict, layer_from_dict, mesh_from_dict
    from .materials import load_material, load_scenarios, scenario_from_dict
    from .simulate import Job, run_job

    with open(a.case) as fh:
        c = yaml.safe_load(fh)
    material = load_material(c["material"])
    scen = (
        [scenario_from_dict(s, material) for s in c["scenarios"]]
        if isinstance(c.get("scenarios"), list)
        else load_scenarios(c.get("scenarios", "configs/scenarios_phase1.yaml"), material)
    )
    job = Job(
        design=design_from_dict(c["design"]),
        seed=int(c.get("seed", 0)),
        layer=layer_from_dict(c.get("layer", {})),
        geometry=geometry_from_dict(c.get("geometry", {})),
        mesh=mesh_from_dict(c.get("mesh", {})),
        material=material,
        scenarios=scen,
        tag="cli",
        contact_model=c.get("contact_model", "series"),
    )
    rows, fails = run_job(job, cache_dir=a.cache)
    for f in fails:
        print(f"FAILED at {f['stage']}: {f['reason']}")
    for r in rows:
        print(f"{r['scenario']:>22s}  R_layer = {r['r_layer']:.4e} m2K/W  k_eff = {r['k_eff']:.3f} W/mK  "
              f"R_total = {r['r_total']:.4e}  balance = {r['heat_balance_error']:.1e}")
    return 0 if rows else 2


def cmd_stage(name):
    def run(a):
        s = _study(a)
        if name == "pilot":
            out = s.pilot()
        elif name == "initial":
            s.initial()
            out = s.surrogate(name="initial")
        elif name == "surrogate":
            import pandas as pd

            out = {"initial": s.surrogate(name="initial")}
            loop = s.out / "search" / "loop_rows.parquet"
            if loop.exists():
                out["loop"] = s.surrogate(pd.read_parquet(loop), name="loop")
        elif name == "search":
            out = s.search()
            import pandas as pd

            s.surrogate(pd.read_parquet(s.out / "search" / "loop_rows.parquet"), name="loop")
        elif name == "final":
            out = s.final()
        elif name == "questions":
            s.q1_q4_grid()
            out = {"q1": s.q1(), "q4": s.q4(), "q2": s.q2(), "q3": s.q3()}
        elif name == "gate":
            out = s.gate()
        else:
            raise ValueError(name)
        print(json.dumps(out, indent=2, default=str)[:4000])
        return 0

    return run


def cmd_figures(a):
    from .figures import make_all

    make_all(_study(a))
    return 0


def cmd_report(a):
    from .report import write_report

    s = _study(a)
    p = write_report(s)
    s.manifest()
    print(f"wrote {p}")
    return 0


def cmd_recommend(a):
    from .recommend import recommend

    res = recommend(
        a.out, material=a.material, thickness=a.thickness, loading=a.loading,
        scenarios=a.scenarios, out_dir=a.export, run_if_unsupported=a.run,
        config=a.config, cache=a.cache, workers=a.workers,
    )
    print(res["text"])
    return 0 if res["supported"] else 3


def cmd_all(a):
    for step in (cmd_verify, cmd_stage("pilot"), cmd_stage("initial"), cmd_stage("search"), cmd_stage("final"),
                 cmd_stage("questions"), cmd_stage("gate"), cmd_figures, cmd_report):
        rc = step(a)
        if rc not in (0, None):
            return rc
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="tcid", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="configs/study_phase1.yaml")
    p.add_argument("--out", default="results/phase1")
    p.add_argument("--cache", default="cache")
    p.add_argument("--workers", type=int, default=None)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("verify").set_defaults(fn=cmd_verify)
    g = sub.add_parser("geometry")
    g.add_argument("--loading", type=float, default=0.3)
    g.add_argument("--d50", default="10 um")
    g.add_argument("--sigma", type=float, default=0.3)
    g.add_argument("--aspect", type=float, default=2.0)
    g.add_argument("--order", type=float, default=0.5)
    g.add_argument("--thickness", default="100 um")
    g.add_argument("--width", default="200 um")
    g.add_argument("--seed", type=int, default=0)
    g.add_argument("--plot", default=None)
    g.set_defaults(fn=cmd_geometry)
    s = sub.add_parser("solve")
    s.add_argument("case")
    s.set_defaults(fn=cmd_solve)
    for name in ("pilot", "initial", "surrogate", "search", "final", "questions", "gate"):
        sub.add_parser(name).set_defaults(fn=cmd_stage(name))
    sub.add_parser("figures").set_defaults(fn=cmd_figures)
    sub.add_parser("report").set_defaults(fn=cmd_report)
    r = sub.add_parser("recommend")
    r.add_argument("--material", default="alumina_epoxy_dev_2d")
    r.add_argument("--thickness", default="100 um")
    r.add_argument("--loading", type=float, default=0.30)
    r.add_argument("--scenarios", default="configs/scenarios_phase1.yaml")
    r.add_argument("--export", default="results/recommendation")
    r.add_argument("--run", action="store_true", help="run a new solver study for unsupported requests")
    r.set_defaults(fn=cmd_recommend)
    sub.add_parser("all").set_defaults(fn=cmd_all)
    a = p.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
