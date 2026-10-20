"""Phase 1 study orchestration (guide sections 5, 6, 8 and 10).

Stages (each reads earlier outputs from disk, so they can be run separately):

    pilot      timing pilot (~20 representative runs) and budget estimate
    initial    space-filling design set with multiple geometry seeds
    surrogate  GP / tree surrogates with design-grouped held-out errors
    search     robust Bayesian optimisation, matched-budget random search, 3-level sweep
    final      fresh seeds + refined meshes for the shortlist and baselines (reported results)
    q1..q4     research questions
    gate       Phase 2 gate assessment
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import replace
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from .config import Design, load_study
from .design_space import DesignSpace
from .materials import Scenario, load_material, load_scenarios, params_for
from .provenance import config_hash, write_manifest
from .search import bayes_opt, best_so_far, design_objectives, job_cost
from .simulate import make_jobs, run_jobs
from .stats import bootstrap, first_crossing, holm, welch
from .surrogate import GPSurrogate, TreeSurrogate, aggregate, cross_validate, features

TARGET_HELDOUT_ERROR = 0.05  # guide section 6 working target


def _json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


class Study:
    def __init__(self, config: str | Path = "configs/study_phase1.yaml", out: str | Path = "results/phase1",
                 cache: str | Path = "cache", workers: int | None = None):
        self.cfg = load_study(config)
        self.root = self.cfg.root
        self.out = Path(out)
        self.cache = Path(cache)
        self.material = load_material(self.cfg.resolve(self.cfg.material_path))
        self.scenarios = load_scenarios(self.cfg.resolve(self.cfg.scenarios_path), self.material)
        self.space = DesignSpace.from_config(self.cfg.design_space)
        self.budget = self.cfg.budget
        self.workers = workers or int(self.budget.get("workers", 4))
        self.base = dict(
            layer=self.cfg.layer,
            geometry=self.cfg.geometry,
            mesh=self.cfg.mesh,
            material=self.material,
            scenarios=self.scenarios,
        )
        self.stress = [s for s in self.scenarios if s.role == "stress"]
        self.nominal = next(s for s in self.scenarios if s.role == "nominal")

    # ------------------------------------------------------------------ io

    def path(self, *parts) -> Path:
        p = self.out.joinpath(*parts)
        p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def save(self, df: pd.DataFrame, *parts):
        df = df.copy()
        df.attrs = {}
        df.to_parquet(self.path(*parts), index=False)

    def load(self, *parts) -> pd.DataFrame:
        return pd.read_parquet(self.out.joinpath(*parts))

    def config_dict(self) -> dict:
        return {
            "study": self.cfg.raw,
            "material": self.material.to_dict(),
            "scenarios": [s.to_dict() for s in self.scenarios],
        }

    # ------------------------------------------------------------------ simulation

    def run(self, designs, seeds, tag, kinds=None, fields=False, **overrides) -> pd.DataFrame:
        base = {**self.base, **overrides}
        meta = (lambda i, d: {"kind": kinds[i]}) if kinds else None
        jobs = make_jobs(designs, seeds, base, tag, meta)
        for j in jobs:
            j.save_fields = fields
        t0 = time.perf_counter()
        rows, fails = run_jobs(jobs, self.workers, self.cache, self.out / "fields" if fields else None)
        print(f"  {tag}: {len(jobs)} jobs, {len(rows)} solves, {len(fails)} failures, {time.perf_counter() - t0:.0f}s", flush=True)
        if len(fails):
            fpath = self.path("failures.parquet")
            prev = pd.read_parquet(fpath) if fpath.exists() else pd.DataFrame()
            allf = pd.concat([prev, fails], ignore_index=True)
            allf = allf.astype({c: str for c in allf.columns if allf[c].dtype == object})
            allf.drop_duplicates(subset=["job_key", "tag"]).to_parquet(fpath, index=False)
        rows = pd.DataFrame(rows)
        failed_ids = set(fails["design_id"]) if len(fails) else set()
        if len(rows):
            # a design counts only if every scenario of every seed solved
            rows = rows[~rows["design_id"].isin(failed_ids)].reset_index(drop=True)
        rows.attrs["failed_design_ids"] = failed_ids
        rows.attrs["n_failures"] = len(fails)
        return rows

    # ------------------------------------------------------------------ stages

    def pilot(self, n: int = 20) -> dict:
        """Time ~20 representative runs including dense and difficult geometries."""
        sp = self.space
        designs = [sp.from_unit(u) for u in sp.sobol(n - 4, 99)]
        b = sp.bounds
        designs += [
            sp.from_unit([0, 1, 1, 0]),  # many small elongated randomly oriented particles
            sp.from_unit([0, 0, 0, 0]),  # many small monodisperse circles
            sp.from_unit([1, 0, 1, 1]),  # large aligned high-aspect particles
            sp.from_unit([1, 1, 0, 0]),  # large polydisperse circles
        ]
        rows = self.run(designs, [0], "pilot")
        self.save(rows, "pilot", "pilot_rows.parquet")
        cost = job_cost(rows)
        fpath = self.out / "failures.parquet"
        nfail = 0
        if fpath.exists():
            f = pd.read_parquet(fpath)
            nfail = int((f["tag"] == "pilot").sum())
        per_job = float(cost["cpu_s"].mean()) if len(cost) else float("nan")
        per_solve = float(rows["solve_time_s"].mean())
        seeds = len(self.budget["seeds"])
        n_init = int(self.budget["initial_designs"])
        n_bo = int(self.budget["bo_iterations"]) * int(self.budget["bo_batch"])
        planned = {
            "initial": n_init * seeds,
            "bayesian_optimisation": n_bo * seeds,
            "random_search": (n_init + n_bo) * seeds,
            "sweep": 3 ** sp.dim * seeds,
        }
        est = {k: v * per_job / self.workers / 3600 for k, v in planned.items()}
        out = {
            "n_runs": int(len(cost)),
            "n_failed": nfail,
            "scenarios_per_job": len(self.scenarios),
            "cpu_s_per_job_mean": per_job,
            "cpu_s_per_job_max": float(cost["cpu_s"].max()),
            "solve_s_mean": per_solve,
            "mesh_s_mean": float(rows.groupby("job_key")["mesh_time_s"].first().mean()),
            "n_elements_mean": float(rows.groupby("job_key")["n_elements"].first().mean()),
            "n_elements_max": float(rows["n_elements"].max()),
            "planned_jobs": planned,
            "planned_solves": {k: v * len(self.scenarios) for k, v in planned.items()},
            "estimated_wall_hours": est,
            "workers": self.workers,
            "decision": "proceed with planned grid" if sum(est.values()) < 6 else "reduce grid",
        }
        _json(self.path("pilot", "pilot.json"), out)
        return out

    def initial(self) -> pd.DataFrame:
        n = int(self.budget["initial_designs"])
        U = self.space.sobol(n, int(self.budget["sobol_seed"]))
        designs = [self.space.from_unit(u) for u in U]
        rows = self.run(designs, self.budget["seeds"], "initial")
        self.save(rows, "design", "initial_rows.parquet")
        return rows

    def surrogate(self, rows: pd.DataFrame | None = None, name: str = "initial") -> dict:
        rows = rows if rows is not None else self.load("design", "initial_rows.parquet")
        agg = aggregate(rows)
        res = {}
        for label, make in (("gaussian_process", lambda: GPSurrogate(n_restarts=2)), ("extra_trees", TreeSurrogate)):
            cv = cross_validate(agg, self.space, make)
            preds = cv.pop("predictions")
            psd = cv.pop("pred_sd")
            res[label] = cv
            if label == "gaussian_process":
                pdf = agg[["design_id", "scenario", "log_r_mean", "log_r_var", "n_seeds"]].copy()
                pdf["pred"] = preds
                pdf["pred_sd"] = psd
                self.save(pdf, "surrogate", f"cv_predictions_{name}.parquet")
        res["meets_target"] = res["gaussian_process"]["median_rel_error"] < TARGET_HELDOUT_ERROR
        res["target"] = TARGET_HELDOUT_ERROR
        _json(self.path("surrogate", f"cv_{name}.json"), res)
        return res

    def search(self) -> dict:
        init = self.load("design", "initial_rows.parquet")
        b = self.budget

        def run(designs, tag, kinds):
            return self.run(designs, b["seeds"], tag, kinds)

        t0 = time.perf_counter()
        bo_rows, hist = bayes_opt(
            run, init, self.space, self.scenarios, self.material,
            int(b["bo_iterations"]), int(b["bo_batch"]), int(b["explore_per_batch"]), int(b["sobol_seed"]),
        )
        t_bo = time.perf_counter() - t0
        self.save(bo_rows, "search", "loop_rows.parquet")
        _json(self.path("search", "bo_history.json"), hist)

        n_total = int(b["initial_designs"]) + int(b["bo_iterations"]) * int(b["bo_batch"])
        U = self.space.uniform(n_total, int(b["random_search_seed"]))
        rnd = self.run([self.space.from_unit(u) for u in U], b["seeds"], "random")
        self.save(rnd, "search", "random_rows.parquet")
        grid = self.space.grid(3)
        sw = self.run([self.space.from_unit(u) for u in grid], b["seeds"], "sweep")
        self.save(sw, "search", "sweep_rows.parquet")
        return self.compare_methods()

    def search_bo_v2(self) -> dict:
        """Corrected-surrogate BO campaign (methods follow-up; plan in docs/gate_v2.md).

        Same initial design, budget, seeds and acquisition settings as the original campaign;
        only the GP noise scaling differs. Its best-of-four is verified and confirmed on the
        existing seed batches and compared with the frozen confirmed picks. The Phase 1
        recommendation is not re-opened.
        """
        b = self.budget
        init = self.load("design", "initial_rows.parquet")

        def run(designs, tag, kinds):
            return self.run(designs, b["seeds"], tag, kinds)

        rows, hist = bayes_opt(
            run, init, self.space, self.scenarios, self.material,
            int(b["bo_iterations"]), int(b["bo_batch"]), int(b["explore_per_batch"]), int(b["sobol_seed"]),
            tag_prefix="bov2",
        )
        self.save(rows, "search", "bo_v2_rows.parquet")
        _json(self.path("search", "bo_v2_history.json"), hist)
        ns = len(self.stress)
        obj = design_objectives(rows, n_stress=ns)
        top = obj.nsmallest(4, "robust")
        designs = {f"bo_v2_{i + 1}": self._design_from_obj(r) for i, (_, r) in enumerate(top.iterrows())}
        ms_ref = replace(self.cfg.mesh, size_factor=float(b["final_refinements"][0]))
        ver = self.run(list(designs.values()), b["final_seeds"], "bov2-verify", mesh=ms_ref)
        stress_names = [s.name for s in self.stress]

        def worst_values(rr, d):
            sel = rr[(rr["design_id"] == config_hash(d.to_dict(), 12)) & rr["scenario"].isin(stress_names)]
            if sel.empty:
                return np.array([])
            worst = sel.groupby("scenario")["r_layer"].mean().idxmax()
            return sel[sel["scenario"] == worst].sort_values("seed")["r_layer"].to_numpy()

        verified = {k: float(worst_values(ver, d).mean()) for k, d in designs.items() if worst_values(ver, d).size}
        pick = min(verified, key=verified.get)
        conf = self.run([designs[pick]], b["confirm_seeds"], "bov2-confirm", mesh=ms_ref)
        pv = worst_values(conf, designs[pick])
        fin = json.loads((self.out / "final" / "final.json").read_text())
        frozen = pd.read_parquet(self.out / "final" / "final_rows.parquet")
        frozen = frozen[frozen["tag"] == "confirm"]
        cands = self._final_designs()
        comps = {}
        for m, lab in fin["method_picks"].items():
            key = lab.split("+")[0]
            other = worst_values(frozen, cands[key])
            comps[f"{m}_minus_bo_v2"] = {**welch(other, pv), "other_label": lab}
        p_adj = holm({k: v["p"] for k, v in comps.items()})
        for k in comps:
            comps[k]["p_holm"] = p_adj[k]
        cv = self.surrogate(rows, name="bo_v2_loop")
        cost = job_cost(rows)
        res = {
            "plan": "docs/gate_v2.md",
            "n_designs": int(obj.shape[0]),
            "solver_calls": int(cost["calls"].sum()),
            "cpu_hours": float(cost["cpu_s"].sum() / 3600),
            "training_acquisition_hours": float(sum(h["train_time_s"] + h["acquisition_time_s"] for h in hist) / 3600),
            "best_robust_loop": float(obj["robust"].min()),
            "top4": {k: d.to_dict() for k, d in designs.items()},
            "verified_robust": verified,
            "pick": pick,
            "confirmed_robust": float(pv.mean()),
            "confirmed_robust_se": float(pv.std(ddof=1) / math.sqrt(pv.size)),
            "comparisons_on_confirmation_batch": comps,
            "cv": {k: v for k, v in cv.items() if k != "target"},
            "note": "methods follow-up; the Phase 1 recommendation is frozen",
        }
        _json(self.path("search", "bo_v2_evaluation.json"), res)
        return res

    def failed_work(self) -> pd.DataFrame:
        """CPU time spent on rejected work per failed job.

        Failure records written since v0.2.0 carry ``elapsed_s``. For older records the rejected
        geometry attempt is re-run (deterministic: same design, seed and settings) and timed.
        """
        fpath = self.out / "failures.parquet"
        if not fpath.exists():
            return pd.DataFrame()
        f = pd.read_parquet(fpath)
        f = f[f["tag"].str.match(r"^(bo-|bov2-|random|sweep)")].drop_duplicates("job_key").copy()
        cache = self.out / "search" / "failed_work.parquet"
        known = pd.read_parquet(cache).set_index("job_key")["elapsed_s"].to_dict() if cache.exists() else {}
        if "elapsed_s" in f:
            known.update(f.dropna(subset=["elapsed_s"]).set_index("job_key")["elapsed_s"].to_dict())
        todo = f[~f["job_key"].isin(known)]
        if len(todo):
            print(f"  re-timing {len(todo)} rejected jobs ...", flush=True)
            args = [(Design(r.loading_target, r.d50, r.sigma_ln, r.aspect, r.order, r.preferred_angle),
                     replace(self.cfg.layer, thickness=r.thickness, width=r.width, margin=r.margin),
                     replace(self.cfg.geometry, min_gap=r.min_gap_setting), int(r.seed)) for r in todo.itertuples()]
            ctx = __import__("multiprocessing").get_context("spawn")
            from concurrent.futures import ProcessPoolExecutor

            with ProcessPoolExecutor(self.workers, mp_context=ctx) as ex:
                times = list(ex.map(_time_generation, args))
            known.update(dict(zip(todo["job_key"], times)))
        f["elapsed_s"] = f["job_key"].map(known)
        pd.DataFrame({"job_key": list(known), "elapsed_s": list(known.values())}).to_parquet(cache, index=False)
        f["method"] = f["tag"].str.replace(r"-\d+$", "", regex=True).map(
            {"bo": "bayesian_optimisation", "bov2": "bo_v2", "random": "random_search", "sweep": "sweep"})
        return f

    def compare_methods(self) -> dict:
        b = self.budget
        loop = self.load("search", "loop_rows.parquet")
        rnd = self.load("search", "random_rows.parquet")
        sw = self.load("search", "sweep_rows.parquet")
        hist = json.loads((self.out / "search" / "bo_history.json").read_text())
        ns = len(self.stress)
        out = {}
        curves = []
        failed = self.failed_work()
        for label, rows in (("bayesian_optimisation", loop), ("random_search", rnd), ("sweep", sw)):
            obj = design_objectives(rows, n_stress=ns)
            cost = job_cost(rows)
            order = list(dict.fromkeys(rows["design_id"]))
            # surrogate/acquisition overhead is charged to BO after the initial design
            overhead = 0.0
            if label == "bayesian_optimisation":
                n_init = rows[rows["tag"] == "initial"]["design_id"].nunique()
                per_it = {f"bo-{h['iteration']}": h["train_time_s"] + h["acquisition_time_s"] for h in hist}
                tags = rows.drop_duplicates("design_id").set_index("design_id")["tag"]
                cur = best_so_far(obj, cost, order)
                extra = np.zeros(len(order))
                seen = set()
                acc = 0.0
                for i, did in enumerate(order):
                    t = tags.get(did, "")
                    if t in per_it and t not in seen:
                        acc += per_it[t]
                        seen.add(t)
                    extra[i] = acc
                cur["cpu_s"] = cur["cpu_s"] + extra
            else:
                cur = best_so_far(obj, cost, order)
            cur["method"] = label
            curves.append(cur)
            n_failed = 0
            fpath = self.out / "failures.parquet"
            if fpath.exists():
                f = pd.read_parquet(fpath)
                prefix = {"bayesian_optimisation": "bo-", "random_search": "random", "sweep": "sweep"}[label]
                n_failed = int(f[f["tag"].str.startswith(prefix)]["design_id"].nunique())
            best = obj.nsmallest(1, "robust").iloc[0]
            fw = failed[failed["method"] == label] if len(failed) else failed
            failed_h = float(fw["elapsed_s"].sum() / 3600) if len(fw) else 0.0
            overhead_h = (sum(h["train_time_s"] + h["acquisition_time_s"] for h in hist) / 3600
                          if label == "bayesian_optimisation" else 0.0)
            out[label] = {
                "n_designs_complete": int(obj.shape[0]),
                "n_designs_requested": int(len(order)) + n_failed,
                "n_designs_rejected": n_failed,
                "solver_calls": int(cost["calls"].sum()),
                "cpu_hours": float(cost["cpu_s"].sum() / 3600),
                "failed_work_cpu_hours": failed_h,
                "surrogate_overhead_hours": overhead_h,
                "total_cpu_hours": float(cost["cpu_s"].sum() / 3600) + failed_h + overhead_h,
                "best_design_id": best["design_id"],
                "best_robust": float(best["robust"]),
                "best_nominal_of_best": float(best["nominal"]),
                "best_design": {k: float(best[k]) for k in ("d50", "sigma_ln", "aspect", "order")},
            }
        curves = pd.concat(curves, ignore_index=True)
        self.save(curves, "search", "best_so_far.parquet")
        overall = min(v["best_robust"] for v in out.values())
        for label in out:
            c = curves[curves["method"] == label]
            hit = c[c["best"] <= overall * 1.01]
            out[label]["calls_to_within_1pct_of_overall_best"] = int(hit["calls"].iloc[0]) if len(hit) else None
            own = out[label]["best_robust"]
            hit_own = c[c["best"] <= own * 1.01]
            out[label]["calls_to_within_1pct_of_own_best"] = int(hit_own["calls"].iloc[0])
        # failure frequency
        fpath = self.out / "failures.parquet"
        if fpath.exists():
            f = pd.read_parquet(fpath)
            fr = f.groupby(f["tag"].str.replace(r"-\d+$", "", regex=True))["stage"].value_counts()
            out["failures_by_tag_and_stage"] = {f"{a}/{b_}": int(v) for (a, b_), v in fr.items()}
        _json(self.path("search", "comparison.json"), out)
        return out

    # ------------------------------------------------------------------ candidates

    def _design_from_obj(self, row) -> Design:
        return Design(
            loading=self.space.loading, d50=float(row["d50"]), sigma_ln=float(row["sigma_ln"]),
            aspect=float(row["aspect"]), order=float(row["order"]), preferred_angle=self.space.preferred_angle,
        )

    # every search method gets the same number of candidates re-verified (no best-of-N advantage)
    METHODS = {
        "bayesian_optimisation": ("robust_opt", "shortlist_2", "shortlist_3", "nominal_opt"),
        "random_search": ("random_search_best", "random_search_2", "random_search_3", "random_search_4"),
        "sweep": ("sweep_best", "sweep_2", "sweep_3", "sweep_4"),
    }
    OPTIMISED = tuple(l for labs in METHODS.values() for l in labs)

    def candidates(self) -> dict[str, Design]:
        """Top four designs of each search method (by loop robust value) and matched baselines.

        The BO four are its three best robust designs plus its best nominal design. Every
        optimised candidate gets its own random-orientation baseline at matched loading,
        sizes and shape (label suffix ``__S0``).
        """
        loop = self.load("search", "loop_rows.parquet")
        ns = len(self.stress)
        obj = design_objectives(loop, n_stress=ns)
        k = int(self.budget.get("shortlist", 3))
        short = obj.nsmallest(k, "robust")
        cands = {}
        for i, (_, r) in enumerate(short.iterrows()):
            cands["robust_opt" if i == 0 else f"shortlist_{i + 1}"] = self._design_from_obj(r)
        nom = obj.nsmallest(1, "nominal").iloc[0]
        cands["nominal_opt"] = self._design_from_obj(nom)
        for method, f in (("random_search", "random_rows.parquet"), ("sweep", "sweep_rows.parquet")):
            o = design_objectives(self.load("search", f), n_stress=ns)
            for lab, (_, r) in zip(self.METHODS[method], o.nsmallest(4, "robust").iterrows()):
                cands[lab] = self._design_from_obj(r)
        for lab in list(cands):
            cands[f"{lab}__S0"] = cands[lab].with_(order=0.0)
        cands["conventional_circles"] = Design(
            loading=self.space.loading, d50=10e-6, sigma_ln=0.3, aspect=1.0, order=0.0,
            preferred_angle=self.space.preferred_angle,
        )
        return cands

    def _summaries(self, rows: pd.DataFrame, labels: dict[str, str]) -> pd.DataFrame:
        ok = rows[rows["status"] == "ok"]
        g = ok.groupby(["design_id", "scenario", "size_factor"]).agg(
            r_mean=("r_layer", "mean"), r_std=("r_layer", "std"), n=("seed", "nunique"),
            k_eff=("k_eff", "mean"), loading=("loading_achieved", "mean"), order_realized=("geo_order_realized", "mean"),
            heat_balance_max=("heat_balance_error", "max"), residual_max=("solver_residual", "max"),
        ).reset_index()
        g["label"] = g["design_id"].map(labels)
        g["r_se"] = g["r_std"] / np.sqrt(g["n"])
        return g

    def final(self) -> dict:
        """Re-simulate the shortlist and baselines with fresh seeds on a refined mesh.

        These solver results, not the loop values or surrogate predictions, are the reported
        results. The recommended design is the verified best optimised candidate.
        """
        from scipy.stats import spearmanr

        from .search import RobustEI
        from .surrogate import physical_features

        b = self.budget
        cands = self.candidates()
        ids: dict[str, list[str]] = {}
        uniq: dict[str, Design] = {}
        for label, d in cands.items():
            did = config_hash(d.to_dict(), 12)
            ids.setdefault(did, []).append(label)
            uniq[did] = d
        labels = {did: "+".join(v) for did, v in ids.items()}
        label_of = {lab: labels[config_hash(d.to_dict(), 12)] for lab, d in cands.items()}
        fresh = b["final_seeds"]
        ms_ref = replace(self.cfg.mesh, size_factor=float(b["final_refinements"][0]))
        rows = self.run(list(uniq.values()), fresh, "final", mesh=ms_ref)
        summ = self._summaries(rows, labels)
        stress_names = [s.name for s in self.stress]

        def stat(did, which):
            s = summ[summ["design_id"] == did]
            if which == "nominal":
                s = s[s["scenario"] == self.nominal.name]
                return (float(s["r_mean"].iloc[0]), float(s["r_se"].iloc[0])) if len(s) else (np.nan, np.nan)
            s = s[s["scenario"].isin(stress_names)]
            if len(s) != len(stress_names):
                return np.nan, np.nan
            w = s.loc[s["r_mean"].idxmax()]
            return float(w["r_mean"]), float(w["r_se"])

        verified = {}
        for did, lab in labels.items():
            r, rse = stat(did, "robust")
            n, nse = stat(did, "nominal")
            verified[lab] = {"design_id": did, "robust": r, "robust_se": rse, "nominal": n, "nominal_se": nse,
                             "design": uniq[did].to_dict()}

        opt_labels = [l for l in self.OPTIMISED if l in cands]
        # selection happens on this verification batch; the confirmation batch below reports
        # values that were not used to choose anything (removes best-of-N selection bias)
        rec = min(opt_labels, key=lambda l: verified[label_of[l]]["robust"])
        rec_d = cands[rec]
        picks = {m: min((l for l in labs if l in cands), key=lambda l: verified[label_of[l]]["robust"])
                 for m, labs in self.METHODS.items()}

        def worst_values(rr, d):
            """Per-realization resistance in the design's worst stress scenario (by mean)."""
            sel = rr[(rr["design_id"] == config_hash(d.to_dict(), 12)) & rr["scenario"].isin(stress_names)]
            worst = sel.groupby("scenario")["r_layer"].mean().idxmax()
            return sel[sel["scenario"] == worst].sort_values("seed")["r_layer"].to_numpy()

        # which optimised candidates did we fail to separate from the recommended one?
        # Welch tests on per-realization worst-scenario values, Holm-corrected over the family.
        rec_vals = worst_values(rows, rec_d)
        tests = {l: welch(worst_values(rows, cands[l]), rec_vals) for l in opt_labels
                 if label_of[l] != label_of[rec]}
        p_adj = holm({l: t["p"] for l, t in tests.items()})
        separable = {l: {**t, "p_holm": p_adj[l], "resolved": bool(p_adj[l] < 0.05 and t["diff"] > 0)}
                     for l, t in tests.items()}

        # confirmation batch: new seeds for the recommended design, every candidate not resolved
        # from it, each method's pick, and the recommended design's matched baseline
        conf_labels = list(dict.fromkeys([rec] + [l for l, v in separable.items() if not v["resolved"]]
                                         + list(picks.values()) + [f"{rec}__S0"]))
        conf_seeds = b["confirm_seeds"]
        conf_designs = {label_of[l]: cands[l] for l in conf_labels}
        conf = self.run(list(conf_designs.values()), conf_seeds, "confirm", mesh=ms_ref)
        conf_summ = self._summaries(conf, labels)
        confirmed = {}
        for lab, d in conf_designs.items():
            vals = worst_values(conf, d)
            nom = conf[(conf["design_id"] == config_hash(d.to_dict(), 12)) & (conf["scenario"] == self.nominal.name)]
            confirmed[lab] = {"robust": float(vals.mean()), "robust_se": float(vals.std(ddof=1) / math.sqrt(vals.size)),
                              "nominal": float(nom["r_layer"].mean()), "n": int(vals.size)}
        rec_conf = worst_values(conf, rec_d)
        conf_tests = {lab: welch(worst_values(conf, d), rec_conf) for lab, d in conf_designs.items()
                      if lab != label_of[rec] and not lab.endswith("__S0")}
        method_vals = {m: worst_values(conf, cands[l]) for m, l in picks.items()}
        method_cmp = {}
        for other in ("random_search", "sweep"):
            t = welch(method_vals[other], method_vals["bayesian_optimisation"])
            base = method_vals[other].mean()
            method_cmp[f"{other}_minus_bayesian_optimisation"] = {
                **t, "rel": t["diff"] / base, "rel_ci_low": t["ci_low"] / base, "rel_ci_high": t["ci_high"] / base}
        rec_vs_base = welch(worst_values(conf, cands[f"{rec}__S0"]), rec_conf)

        # paired candidate - matched-baseline differences under the same scenarios
        diffs = []
        for l in opt_labels:
            did, bid = config_hash(cands[l].to_dict(), 12), config_hash(cands[f"{l}__S0"].to_dict(), 12)
            for sc in summ["scenario"].unique():
                cv = rows[(rows["design_id"] == did) & (rows["scenario"] == sc)]["r_layer"].to_numpy()
                bv = rows[(rows["design_id"] == bid) & (rows["scenario"] == sc)]["r_layer"].to_numpy()
                if not cv.size or not bv.size:
                    continue
                t = welch(bv, cv)
                diffs.append({"label": l, "scenario": sc, "baseline_minus_candidate": t["diff"], "se": t["se"],
                              "ci_low": t["ci_low"], "ci_high": t["ci_high"], "p": t["p"],
                              "relative_benefit": t["diff"] / bv.mean(), "t_ratio": t["diff"] / t["se"] if t["se"] else float("inf")})
        diffs = pd.DataFrame(diffs)
        self.save(diffs, "final", "paired_differences.parquet")

        # fields and mesh refinement for the key claim: recommended design vs its baseline
        key = [rec_d, cands[f"{rec}__S0"]]
        roles = {config_hash(rec_d.to_dict(), 12): "recommended", config_hash(key[1].to_dict(), 12): "baseline"}
        fld = self.run(key, fresh[:1], "final-fields", mesh=ms_ref, fields=True, scenarios=[self.nominal])
        ms_fine = replace(self.cfg.mesh, size_factor=float(b["final_refinements"][1]))
        fine = self.run(key, fresh[:2], "final-refine", mesh=ms_fine)
        coarse = self.run(key, fresh[:2], "final-coarse", mesh=self.cfg.mesh)
        allrows = pd.concat([rows, fine, coarse, fld, conf], ignore_index=True)
        self.save(allrows, "final", "final_rows.parquet")
        self.save(summ, "final", "final_summary.parquet")
        self.save(conf_summ, "final", "confirmation_summary.parquet")
        ref = []
        for did, role in roles.items():
            for sc in [s.name for s in self.scenarios]:
                vals = {}
                for tag, rr in (("loop", coarse), ("refined", rows), ("fine", fine)):
                    sel = rr[(rr["design_id"] == did) & (rr["scenario"] == sc) & rr["seed"].isin(fresh[:2])]
                    vals[tag] = float(sel["r_layer"].mean()) if len(sel) else float("nan")
                ref.append({"design_id": did, "role": role, "label": labels[did], "scenario": sc, **vals,
                            "change_loop_to_refined": abs(vals["refined"] - vals["loop"]) / vals["refined"],
                            "change_refined_to_fine": abs(vals["fine"] - vals["refined"]) / vals["fine"]})
        ref = pd.DataFrame(ref)
        self.save(ref, "final", "mesh_refinement.parquet")
        mesh_err = float(ref["change_refined_to_fine"].max())

        out_c = {lab: d.to_dict() for lab, d in cands.items()}
        out_c["recommended"] = rec_d.to_dict()
        out_c["baseline_random_orientation"] = key[1].to_dict()
        _json(self.path("final", "candidates.json"), out_c)

        # surrogate trained on all loop data: does its ordering survive direct checks?
        loop = self.load("search", "loop_rows.parquet")
        agg = aggregate(loop)
        gp = GPSurrogate(n_restarts=3).fit(features(agg, self.space), agg["log_r_mean"].to_numpy(),
                                            (agg["log_r_var"] / agg["n_seeds"].clip(lower=1)).to_numpy())
        joblib.dump({"model": gp.model, "space": self.space,
                     "features": ["d50 (log-scaled)", "sigma_ln", "aspect", "order", "log10 R_interface (scaled)",
                                  "k_filler (scaled)"],
                     "target": "log R_layer [m2K/W]", "trained_on": "initial + BO designs, loop mesh, training seeds"},
                    self.path("surrogate", "gp_final.joblib"))
        phys = physical_features(np.array([s.r_interface for s in self.stress]),
                                 np.array([self.material.k_filler_par * s.k_filler_scale for s in self.stress]))
        acq = RobustEI(gp, self.space, phys)
        pred = {}
        for did, lab in labels.items():
            mu, _ = acq.predict(self.space.to_unit(uniq[did])[None, :])
            pred[lab] = float(np.exp(mu).max())
        order_pred = sorted(pred, key=pred.get)
        order_ver = sorted(verified, key=lambda k: verified[k]["robust"])
        rho = spearmanr([pred[k] for k in order_ver], [verified[k]["robust"] for k in order_ver]).statistic
        obj = design_objectives(loop, n_stress=len(self.stress))
        rho_nr = spearmanr(obj["nominal"], obj["robust"]).statistic
        loop_top = verified[label_of["robust_opt"]]
        rec_benefit = diffs[diffs["label"] == rec]["relative_benefit"]
        res = {
            "seeds": fresh,
            "mesh_size_factor": ms_ref.size_factor,
            "verified": verified,
            "recommended": rec,
            "recommended_label": label_of[rec],
            "loop_top1": "robust_opt",
            "loop_top1_value_on_training_seeds": float(obj["robust"].min()),
            "loop_top1_value_on_fresh_seeds": loop_top["robust"],
            "selection_batch_seeds": fresh,
            "confirmation_seeds": conf_seeds,
            "separability_from_recommended": separable,
            "method_picks": {m: label_of[l] for m, l in picks.items()},
            "confirmed": confirmed,
            "confirmed_vs_recommended": conf_tests,
            "confirmed_recommended_vs_baseline": {**rec_vs_base, "rel": rec_vs_base["diff"] / worst_values(conf, cands[f"{rec}__S0"]).mean()},
            "method_comparison_confirmed": method_cmp,
            "surrogate_predicted_robust": pred,
            "order_predicted": order_pred,
            "order_verified": order_ver,
            "top1_survives": label_of[rec] == label_of["robust_opt"],
            "top1_not_resolved_from_verified_best": (label_of["robust_opt"] == label_of[rec])
            or not separable["robust_opt"]["resolved"],
            "spearman_pred_vs_verified": float(rho),
            "spearman_nominal_vs_robust_loop": float(rho_nr),
            "mesh_refinement_max_change": mesh_err,
            "mesh_refinement_loop_to_refined_max": float(ref["change_loop_to_refined"].max()),
            "smallest_relative_benefit_claimed": float(rec_benefit.min()) if len(rec_benefit) else None,
            "numerical_error_below_claim": bool(len(rec_benefit) == 0 or mesh_err < 0.2 * rec_benefit.min()),
            "labels": labels,
            "label_of": label_of,
            "nominal_vs_robust_same_design": label_of["robust_opt"] == label_of["nominal_opt"],
        }
        _json(self.path("final", "final.json"), res)
        return res

    # ------------------------------------------------------------------ research questions

    def _final_designs(self) -> dict[str, Design]:
        c = json.loads((self.out / "final" / "candidates.json").read_text())
        return {k: Design(**v) for k, v in c.items()}

    def _scenario_grid(self, r_values, kf_values) -> list[Scenario]:
        out = []
        for r in r_values:
            for kf in kf_values:
                role = "boundary" if (r > 1e-6 or r < 1e-8) else "stress"
                out.append(Scenario(name=f"R{r:.0e}_kf{kf:g}", r_interface=r, k_filler_scale=kf / self.material.k_filler_par, role=role))
        return out

    def q1_q4_grid(self) -> pd.DataFrame:
        """Shared grid for Q1 (design map) and Q4 (value of information)."""
        c = self._final_designs()
        ro = c["recommended"]
        lo, hi = self.space.bounds["d50"]
        # size ladder sharing the recommended width, shape and alignment, plus the recommended design itself
        cands = {
            "baseline_random_orientation": c["baseline_random_orientation"],
            "recommended": ro,
            **self._size_ladder(ro),
        }
        r_values = [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5]
        scen = self._scenario_grid(r_values, [20.0, 30.0, 35.0])
        designs, meta = [], []
        for lab, d in cands.items():
            fracs = [1.0] if d.order == 0 else [1.0, 0.75, 0.5]
            for f in fracs:
                designs.append(d.with_(order=d.order * f))
                meta.append({"label": lab, "order_fraction": f})
        seeds = [2000, 2001, 2002]
        rows = self.run(designs, seeds, "q1q4", scenarios=scen)
        did_meta = {config_hash(d.to_dict(), 12): m for d, m in zip(designs, meta)}
        rows["label"] = rows["design_id"].map(lambda x: did_meta[x]["label"])
        rows["order_fraction"] = rows["design_id"].map(lambda x: did_meta[x]["order_fraction"])
        self.save(rows, "questions", "q1q4_rows.parquet")
        return rows

    def q1(self) -> dict:
        rows = self.load("questions", "q1q4_rows.parquet")
        ok = rows[rows["status"] == "ok"]
        g = ok.groupby(["label", "order_fraction", "r_interface", "k_filler_par"]).agg(
            r_mean=("r_layer", "mean"), r_std=("r_layer", "std"), n=("seed", "nunique")).reset_index()
        g["se"] = g["r_std"] / np.sqrt(g["n"])
        base = g[g["label"] == "baseline_random_orientation"].set_index(["r_interface", "k_filler_par"])

        def vals(label, f, r, kf):
            f = 1.0 if label == "baseline_random_orientation" else f
            sel = ok[(ok["label"] == label) & (ok["order_fraction"] == f) & (ok["r_interface"] == r) & (ok["k_filler_par"] == kf)]
            return sel["r_layer"].to_numpy()

        recs = []
        for (f, r, kf), sub in g[g["label"] != "baseline_random_orientation"].groupby(["order_fraction", "r_interface", "k_filler_par"]):
            bb = base.loc[(r, kf)]
            pool = pd.concat([sub, bb.to_frame().T.assign(order_fraction=f, r_interface=r, k_filler_par=kf)])
            pool = pool.sort_values("r_mean")
            w, ru = pool.iloc[0], pool.iloc[1]
            t = welch(vals(ru["label"], f, r, kf), vals(w["label"], f, r, kf))
            recs.append({"order_fraction": f, "r_interface": r, "k_filler": kf, "winner": w["label"], "runner_up": ru["label"],
                         "margin": t["diff"], "margin_rel": t["diff"] / w["r_mean"], "p_welch": t["p"],
                         "distinguishable": bool(t["p"] < 0.05),
                         "baseline_rank": int(list(pool["label"]).index("baseline_random_orientation")) + 1})
        m = pd.DataFrame(recs)
        self.save(g, "questions", "q1_summary.parquet")
        self.save(m, "questions", "q1_map.parquet")
        supported = m[(m["r_interface"] >= 1e-8) & (m["r_interface"] <= 1e-6)]
        res = {
            "winners_supported_range": supported["winner"].value_counts().to_dict(),
            "winners_all": m["winner"].value_counts().to_dict(),
            "baseline_ever_wins_supported": bool((supported["winner"] == "baseline_random_orientation").any()),
            "indistinguishable_fraction_supported": float(1 - supported["distinguishable"].mean()),
            "seeds_per_cell": int(g["n"].min()),
            "cell_test": "Welch two-sample test of winner vs runner-up, p < 0.05 (screening only, not Holm-corrected)",
            "coarse_grid_crossings": self._reversals(g),
            "note": "coarse-grid crossings are crossings of sample means at 3 seeds; see q1_crossings.json for "
                    "dense-grid, 12-seed bootstrap intervals",
        }
        _json(self.path("questions", "q1.json"), res)
        return res

    @staticmethod
    def _reversals(g: pd.DataFrame) -> list[dict]:
        """Interface resistance at which the ordering of two candidates changes (k_f = 30, full order)."""
        out = []
        sub = g[(g["k_filler_par"] == 30.0) & (g["order_fraction"] == 1.0)]
        labs = sorted(sub["label"].unique())
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                a = sub[sub["label"] == labs[i]].set_index("r_interface")["r_mean"]
                b = sub[sub["label"] == labs[j]].set_index("r_interface")["r_mean"]
                d = (a - b).dropna().sort_index()
                s = np.sign(d.to_numpy())
                for k in range(len(s) - 1):
                    if s[k] != s[k + 1] and s[k] != 0:
                        x0, x1 = math.log10(d.index[k]), math.log10(d.index[k + 1])
                        y0, y1 = d.iloc[k], d.iloc[k + 1]
                        cross = 10 ** (x0 + (x1 - x0) * y0 / (y0 - y1))  # log-linear interpolation
                        out.append({"a": labs[i], "b": labs[j], "between": [float(d.index[k]), float(d.index[k + 1])],
                                    "crossing_estimate": float(cross),
                                    "better_below": labs[j] if y0 > 0 else labs[i]})
        return out

    def feasibility(self) -> dict:
        """Separate generator limitations from geometric impossibility for the excluded cases.

        For each case excluded after random sequential placement failed, the deterministic shelf
        construction (geometry.construct_shelf) is tried with the same population and
        orientations. Success proves that a valid arrangement exists; it is an ordered
        arrangement, so it is used only as a feasibility statement plus one labelled
        illustrative solve, never mixed with random realizations.
        """
        from .geometry import FAILURE_KINDS, construct_shelf, generate, validate
        from .mesh import mesh_particles
        from .solver import LayerModel

        c = self._final_designs()
        ro, med = c["recommended"], self._size_ladder(c["recommended"])["medium_particles_aligned"]
        t0 = self.cfg.layer.thickness
        thin = replace(self.cfg.layer, thickness=0.5 * t0, width=t0)
        cases = {
            "small rung 4.0 um (recommended shape)": (ro.with_(d50=4.0e-6), self.cfg.layer, self.cfg.geometry, [2000, 2001, 2002]),
            "small rung 4.5 um (recommended shape)": (ro.with_(d50=4.5e-6), self.cfg.layer, self.cfg.geometry, [2000, 2001, 2002]),
            "recommended at 50 um": (ro, thin, self.cfg.geometry, self.budget["confirm_seeds"][:3]),
            "medium pair candidate at 50 um": (med, thin, self.cfg.geometry, self.budget["confirm_seeds"]),
            "medium pair baseline at 50 um": (med.with_(order=0.0), thin, self.cfg.geometry, self.budget["confirm_seeds"]),
            "recommended with 6 um margin": (ro, replace(self.cfg.layer, margin=6e-6), self.cfg.geometry,
                                             self.budget.get("perturbation_seeds", [4000])),
        }
        out = {"failure_kinds": FAILURE_KINDS, "cases": {}}
        for name, (d, lay, gs, seeds) in cases.items():
            rec = {"rsa": {}, "shelf": {}}
            for sd in seeds:
                g = generate(d, lay, gs, sd)
                rec["rsa"][str(sd)] = "ok" if g.ok else g.reason
                h = construct_shelf(d, lay, gs, sd)
                v = validate(h) if h.ok else {"valid": False}
                rec["shelf"][str(sd)] = "valid" if (h.ok and v["valid"]) else (h.reason or "invalid")
            n_rsa = sum(v == "ok" for v in rec["rsa"].values())
            n_shelf = sum(v == "valid" for v in rec["shelf"].values())
            geometric = all("cannot fit" in v or "truncated" in v for v in rec["shelf"].values()) and n_rsa == 0
            if n_rsa == len(seeds):
                verdict = "feasible with random placement for every realization"
            elif n_rsa > 0 or n_shelf > 0:
                verdict = (f"feasible: valid arrangements exist ({n_rsa}/{len(seeds)} by random placement, "
                           f"{n_shelf}/{len(seeds)} by shelf construction); the remaining failures are generator "
                           "limitations, not shown impossible")
            elif geometric:
                verdict = "impossible under the stated size rules (particles longer than the layer allows)"
            else:
                verdict = "not shown feasible by either strategy (neither proves impossibility)"
            rec.update({"n_seeds": len(seeds), "rsa_ok": n_rsa, "shelf_ok": n_shelf, "verdict": verdict})
            out["cases"][name] = rec
        # illustrative solve of a constructed arrangement for the smallest excluded size
        d = ro.with_(d50=4.0e-6)
        ill = []
        for sd in (2000, 2001, 2002):
            h = construct_shelf(d, self.cfg.layer, self.cfg.geometry, sd)
            if not (h.ok and validate(h)["valid"]):
                continue
            lm = LayerModel(mesh_particles(h, self.cfg.mesh))
            for sc in (self.nominal, max(self.stress, key=lambda x: (x.r_interface, -x.k_filler_scale))):
                ill.append({"seed": sd, "scenario": sc.name, "r_layer": lm.solve(params_for(self.material, sc)).r_layer})
        out["illustrative_constructed_4um"] = {
            "note": "ordered shelf arrangement; NOT comparable with randomly placed realizations",
            "solves": ill}
        # does conditioning on successful placement distort the realized orientation distribution?
        acc = pd.concat([self.load("search", f) for f in ("loop_rows.parquet", "random_rows.parquet", "sweep_rows.parquet")])
        acc = acc.drop_duplicates("job_key")
        acc = acc[acc["aspect"] > 1.2]
        dd = acc["geo_order_realized"] - acc["order"]
        out["orientation_bias_check"] = {
            "accepted_geometries": int(len(acc)),
            "mean_realized_minus_target_order": float(dd.mean()),
            "se": float(dd.std(ddof=1) / math.sqrt(len(dd))),
            "geometries_needing_restarts": int((acc["geo_placement_restarts"] > 0).sum()),
        }
        fpath = self.out / "failures.parquet"
        if fpath.exists():
            from .geometry import classify_failure

            f = pd.read_parquet(fpath)
            f["failure_kind"] = [classify_failure(a, b) for a, b in zip(f["stage"], f["reason"])]
            out["failure_counts_by_kind"] = f.groupby("failure_kind")["job_key"].nunique().to_dict()
        _json(self.path("questions", "feasibility.json"), out)
        return out

    def _size_ladder(self, ro: Design) -> dict[str, Design]:
        """Small / medium / large particles sharing the recommended width, shape and alignment.

        The small rung is the smallest size (from the lower bound, in 0.5 um steps) that packs at
        the study loading for three test seeds; infeasible smaller sizes are recorded, not hidden.
        """
        from .geometry import generate

        lo, hi = self.space.bounds["d50"]
        small, rejected = lo, []
        while small < hi:
            if all(generate(ro.with_(d50=small), self.cfg.layer, self.cfg.geometry, sd).ok for sd in (2000, 2001, 2002)):
                break
            rejected.append(small)
            small += 0.5e-6
        self._ladder_rejected = rejected
        return {"small_particles_aligned": ro.with_(d50=small),
                "medium_particles_aligned": ro.with_(d50=math.sqrt(small * hi)),
                "large_particles_aligned": ro.with_(d50=hi * 0.85)}

    def q1_crossings(self) -> dict:
        """Targeted runs around the particle-size crossings with bootstrap intervals on their location.

        Dense interface-resistance grid (8 points per decade) at k_filler = 30 W/mK and full
        orientation; 12 seeds per design. Each bootstrap resample redraws seeds independently per
        design and locates the first sign change of the mean-difference curve.
        """
        c = self._final_designs()
        designs = self._size_ladder(c["recommended"])
        r_grid = list(np.logspace(-7, -5, 17))
        scen = self._scenario_grid(r_grid, [30.0])
        seeds = self.budget["crossing_seeds"]
        rows = self.run(list(designs.values()), seeds, "q1-crossings", scenarios=scen)
        lab = {config_hash(d.to_dict(), 12): k for k, d in designs.items()}
        rows["label"] = rows["design_id"].map(lab)
        self.save(rows, "questions", "q1_crossing_rows.parquet")
        ok = rows[rows["status"] == "ok"]
        curves = {}
        for k in designs:
            p = ok[ok["label"] == k].pivot_table(index="seed", columns="r_interface", values="r_layer")
            curves[k] = p.dropna()
        x = np.array(sorted(curves["small_particles_aligned"].columns))
        out = {"seeds": seeds, "r_grid": [float(v) for v in x],
               "d50_um": {k: d.d50 * 1e6 for k, d in designs.items()},
               "infeasible_small_sizes_um": [v * 1e6 for v in getattr(self, "_ladder_rejected", [])], "pairs": {}}
        rec_rows = []
        for a, b in (("small_particles_aligned", "medium_particles_aligned"), ("medium_particles_aligned", "large_particles_aligned"),
                     ("small_particles_aligned", "large_particles_aligned")):
            if len(curves[a]) < 3 or len(curves[b]) < 3:
                out["pairs"][f"{a}|{b}"] = {"status": "not enough feasible realizations"}
                continue
            A, B = curves[a][x].to_numpy(), curves[b][x].to_numpy()
            bs = bootstrap({"A": A, "B": B}, lambda A, B: first_crossing(x, A - B), n=4000, seed=3)
            dmean = A.mean(0) - B.mean(0)
            supported_hi = 1e-6
            if not np.isfinite(bs["estimate"]):
                status = "no crossing on the grid"
            elif bs["ci_high"] <= supported_hi:
                status = "inside supported range"
            elif bs["ci_low"] > supported_hi:
                status = "above supported range (boundary stress case)"
            else:
                status = "interval straddles the supported-range limit (1e-6)"
            out["pairs"][f"{a}|{b}"] = {**bs, "better_below": a if dmean[0] < 0 else b, "status": status,
                                        "per_seed_mean_difference_at_grid": [float(v) for v in dmean]}
            for i, xv in enumerate(x):
                dd = bootstrap({"A": A[:, i], "B": B[:, i]}, lambda A, B: A - B, n=1000, seed=4)
                rec_rows.append({"pair": f"{a}|{b}", "r_interface": xv, "diff": dd["estimate"], "ci_low": dd["ci_low"],
                                 "ci_high": dd["ci_high"]})
        self.save(pd.DataFrame(rec_rows), "questions", "q1_crossing_curves.parquet")
        _json(self.path("questions", "q1_crossings.json"), out)
        return out

    def q4(self) -> dict:
        """Sensitivity shares: how much of the spread of a decision quantity disappears when one
        input is known exactly, averaged over a uniform grid of the *chosen* input ranges.

        This is a sensitivity index conditional on the ranges, grid weighting and independence
        assumptions; it is not a value-of-information analysis (no measurement cost, accuracy
        or decision change is modelled). Range variants show how the shares move.
        """
        rows = self.load("questions", "q1q4_rows.parquet")
        ok = rows[rows["status"] == "ok"]
        g_all = ok.groupby(["label", "order_fraction", "r_interface", "k_filler_par"])["r_layer"].mean().reset_index()
        inputs = {"orientation_achieved": "order_fraction", "interface_resistance": "r_interface",
                  "filler_conductivity": "k_filler_par"}

        def shares(r_lo, r_hi, fracs):
            g = g_all[(g_all["r_interface"] >= r_lo) & (g_all["r_interface"] <= r_hi)]
            base = g[g["label"] == "baseline_random_orientation"].set_index(["r_interface", "k_filler_par"])["r_layer"]
            out = {}
            for lab in ("recommended", "small_particles_aligned", "medium_particles_aligned", "large_particles_aligned"):
                c = g[(g["label"] == lab) & g["order_fraction"].isin(fracs)].copy()
                if c.empty:
                    continue
                c["delta"] = c.apply(lambda r: base.loc[(r["r_interface"], r["k_filler_par"])] - r["r_layer"], axis=1)
                out[lab] = _variance_reduction(c, "delta", inputs)
            a = g[(g["label"] == "small_particles_aligned") & g["order_fraction"].isin(fracs)].set_index(
                ["order_fraction", "r_interface", "k_filler_par"])["r_layer"]
            bq = g[(g["label"] == "large_particles_aligned") & g["order_fraction"].isin(fracs)].set_index(
                ["order_fraction", "r_interface", "k_filler_par"])["r_layer"]
            dd = (bq - a).dropna().reset_index(name="delta")
            out["large_minus_small"] = _variance_reduction(dd, "delta", inputs)
            return out

        variants = {
            "reference: R 1e-8..1e-6, achieved/target order 0.5-1.0": (1e-8, 1e-6, [0.5, 0.75, 1.0]),
            "narrower orientation uncertainty (0.75-1.0)": (1e-8, 1e-6, [0.75, 1.0]),
            "wider interface range (1e-8..1e-5, includes stress cases)": (1e-8, 1e-5, [0.5, 0.75, 1.0]),
            "narrower interface range (1e-8..1e-7)": (1e-8, 1e-7, [0.5, 0.75, 1.0]),
        }
        res = {"definition": "1 - E[Var(delta | input)] / Var(delta) over a uniform grid of the stated ranges; "
                             "grid means over 3 seeds (geometry noise included in the spread)",
               "variants": {k: shares(*v) for k, v in variants.items()}}
        ref = next(iter(res["variants"].values()))
        res.update(ref)  # reference shares at top level (backwards compatible)
        _json(self.path("questions", "q4.json"), res)
        return res

    def q2(self) -> dict:
        """Contact-resistance decision boundary with bootstrap intervals.

        Uses the confirmation seeds and the refined (reported) mesh, so the reference-thickness
        values of the recommended pair are exactly the confirmation-batch values in the report.
        A second pair (medium-size ladder design vs its own random-orientation baseline) is run
        at every thickness so the thickness trend can be tested even where the recommended
        distribution does not fit in the layer.
        """
        c = self._final_designs()
        med = self._size_ladder(c["recommended"])["medium_particles_aligned"]
        pairs = {"recommended": (c["recommended"], c["baseline_random_orientation"]),
                 "medium_particles_aligned": (med, med.with_(order=0.0))}
        seeds = self.budget["confirm_seeds"]
        ms_ref = replace(self.cfg.mesh, size_factor=float(self.budget["final_refinements"][0]))
        t0 = self.cfg.layer.thickness
        lab = {}
        for pname, (cd, bd) in pairs.items():
            lab[config_hash(cd.to_dict(), 12)] = (pname, "candidate")
            lab[config_hash(bd.to_dict(), 12)] = (pname, "baseline")
        designs = [d for pr in pairs.values() for d in pr]
        rows = []
        for f in (0.5, 1.0, 1.5):
            layer = self.cfg.layer if f == 1.0 else replace(self.cfg.layer, thickness=f * t0, width=2 * f * t0)
            rows.append(self.run(designs, seeds, f"q2-t{f * t0 * 1e6:.0f}", layer=layer, mesh=ms_ref))
        rows = pd.concat(rows, ignore_index=True)
        rows["pair"] = rows["design_id"].map(lambda x: lab[x][0])
        rows["role"] = rows["design_id"].map(lambda x: lab[x][1])
        # Robin check of the series reduction at the nominal interface (same seeds and mesh)
        robin_sc = [replace(self.nominal, name=f"robin_{rc:.0e}", r_top=rc, r_bottom=rc) for rc in (1e-5, 1e-4)]
        rb = self.run(designs, seeds[:2], "q2-robin", scenarios=robin_sc, contact_model="robin", mesh=ms_ref)
        rb["pair"] = rb["design_id"].map(lambda x: lab[x][0])
        rb["role"] = rb["design_id"].map(lambda x: lab[x][1])
        self.save(rows, "questions", "q2_rows.parquet")
        self.save(rb, "questions", "q2_robin_rows.parquet")
        ok = rows[rows["status"] == "ok"]
        contacts = np.array([0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4])
        recs, bounds, layer_tab, zero = [], {}, [], {}
        infeasible = {}
        for pname in pairs:
            for f in (0.5, 1.0, 1.5):
                t = f * t0
                tkey = f"{t * 1e6:.0f}um"
                sub_t = ok[(ok["pair"] == pname) & np.isclose(ok["thickness"], t)]
                if sub_t["role"].nunique() < 2:
                    infeasible.setdefault(pname, []).append(tkey)
                    continue
                for sc, sub in sub_t.groupby("scenario"):
                    bv = sub[sub["role"] == "baseline"].sort_values("seed")["r_layer"].to_numpy()
                    cv = sub[sub["role"] == "candidate"].sort_values("seed")["r_layer"].to_numpy()
                    if bv.size < 3 or cv.size < 3:
                        continue
                    layer_tab.append({"pair": pname, "thickness": t, "scenario": sc, "r_baseline": bv.mean(),
                                      "r_baseline_se": bv.std(ddof=1) / math.sqrt(bv.size), "r_candidate": cv.mean(),
                                      "r_candidate_se": cv.std(ddof=1) / math.sqrt(cv.size),
                                      "n_baseline": bv.size, "n_candidate": cv.size})
                    for rct in contacts:
                        bs = bootstrap({"b": bv, "c": cv}, lambda b, c, rct=rct: (b - c) / (b + rct), n=2000, seed=1)
                        recs.append({"pair": pname, "thickness": t, "scenario": sc, "r_contact_total": rct,
                                     "r_layer_baseline": bv.mean(), "r_layer_candidate": cv.mean(),
                                     "benefit_abs": bv.mean() - cv.mean(), "benefit_rel_total": bs["estimate"],
                                     "benefit_ci_low": bs["ci_low"], "benefit_ci_high": bs["ci_high"],
                                     "contact_share_baseline": rct / (bv.mean() + rct)})
                        if rct == 0:
                            zero.setdefault(pname, {})[f"{tkey}/{sc}"] = bs
                    for thr in (0.10, 0.05):
                        bs = bootstrap({"b": bv, "c": cv}, lambda b, c, thr=thr: max(0.0, (b - c) / thr - b), n=4000, seed=2)
                        bounds.setdefault(f"{thr:.2f}", {}).setdefault(pname, {}).setdefault(tkey, {})[sc] = bs
        q2 = pd.DataFrame(recs)
        self.save(q2, "questions", "q2_benefit.parquet")
        self.save(pd.DataFrame(layer_tab), "questions", "q2_layer_resistances.parquet")
        # series reduction check: R_top + R_layer(isothermal faces) + R_bottom vs the Robin assembly solve
        rbok = rb[rb["status"] == "ok"]
        iso = ok[np.isclose(ok["thickness"], t0) & (ok["scenario"] == self.nominal.name)].set_index(["design_id", "seed"])["r_layer"]
        recs_r = []
        for _, r in rbok.iterrows():
            r_iso = iso.get((r["design_id"], r["seed"]))
            if r_iso is None:
                continue
            series = r["r_top"] + r_iso + r["r_bottom"]
            recs_r.append({"pair": r["pair"], "role": r["role"], "scenario": r["scenario"], "seed": r["seed"],
                           "r_contact_each": r["r_top"], "r_total_robin": r["r_total"], "r_total_series": series,
                           "series_error_rel": (series - r["r_total"]) / r["r_total"]})
        robin = pd.DataFrame(recs_r)
        self.save(robin, "questions", "q2_robin.parquet")
        # realizations rejected per pair, role and thickness (a design with any rejection is excluded)
        rejected = {}
        fpath = self.out / "failures.parquet"
        if fpath.exists():
            fl = pd.read_parquet(fpath)
            fl = fl[fl["tag"].str.startswith("q2-t")]
            for did, (pname, role) in lab.items():
                for tag, n in fl[fl["design_id"] == did].groupby("tag")["seed"].nunique().items():
                    rejected.setdefault(pname, {}).setdefault(role, {})[tag.replace("q2-t", "") + "um"] = (
                        f"{n}/{len(seeds)} realizations rejected: "
                        + str(fl[(fl["design_id"] == did) & (fl["tag"] == tag)]["reason"].iloc[0])[:90])
        res = {
            "seeds": seeds,
            "mesh_size_factor": ms_ref.size_factor,
            "pairs": {k: {"candidate": v[0].to_dict(), "baseline": v[1].to_dict()} for k, v in pairs.items()},
            "rejected_realizations": rejected,
            "infeasible_thicknesses": infeasible,
            "contact_total_where_benefit_falls_to": bounds,
            "benefit_at_zero_contact": zero,
            "robin_vs_series_max_rel_diff": float(robin["series_error_rel"].abs().max()) if len(robin) else None,
            "robin_vs_series_by_pair": robin.groupby("pair")["series_error_rel"].apply(lambda x: float(x.abs().max())).to_dict()
            if len(robin) else {},
            "note": "Identical series contacts cannot reorder candidates; they reduce the relative benefit. "
                    "Thresholds amplify noise ~10x (threshold = dR/0.1 - R_baseline), hence the bootstrap intervals.",
        }
        _json(self.path("questions", "q2.json"), res)
        return res

    def q3(self) -> dict:
        c = self._final_designs()
        ro, base = c["recommended"], c["baseline_random_orientation"]
        lay, gs = self.cfg.layer, self.cfg.geometry
        variants = {
            "reference": ({}, ro, base),
            "order_minus_0.1": ({}, ro.with_(order=max(0.0, ro.order - 0.1)), None),
            "order_minus_0.2": ({}, ro.with_(order=max(0.0, ro.order - 0.2)), None),
            "d50_x0.8": ({}, ro.with_(d50=ro.d50 * 0.8), base.with_(d50=base.d50 * 0.8)),
            "d50_x1.2": ({}, ro.with_(d50=ro.d50 * 1.2), base.with_(d50=base.d50 * 1.2)),
            "sigma_plus_0.1": ({}, ro.with_(sigma_ln=ro.sigma_ln + 0.1), base.with_(sigma_ln=base.sigma_ln + 0.1)),
            "min_gap_0.5um": ({"geometry": replace(gs, min_gap=0.5e-6)}, ro, base),
            "min_gap_2um": ({"geometry": replace(gs, min_gap=2.0e-6)}, ro, base),
            "margin_1.5um": ({"layer": replace(lay, margin=1.5e-6)}, ro, base),
            "margin_6um": ({"layer": replace(lay, margin=6e-6)}, ro, base),
            "width_x2": ({"layer": replace(lay, width=2 * lay.width)}, ro, base),
        }
        seeds = self.budget.get("perturbation_seeds", [4000, 4001, 4002, 4003])
        allrows = []
        for name, (ov, dc, db) in variants.items():
            ds = [dc] + ([db] if db is not None else [])
            r = self.run(ds, seeds, f"q3-{name}", **ov)
            r["variant"] = name
            r["role"] = r["design_id"].map({config_hash(dc.to_dict(), 12): "candidate",
                                            **({config_hash(db.to_dict(), 12): "baseline"} if db is not None else {})})
            allrows.append(r)
        rows = pd.concat(allrows, ignore_index=True)
        self.save(rows, "questions", "q3_rows.parquet")
        ok = rows[rows["status"] == "ok"]
        g = ok.groupby(["variant", "role", "scenario"]).agg(r=("r_layer", "mean"), sd=("r_layer", "std"),
                                                            n=("seed", "nunique")).reset_index()
        g["se"] = g["sd"] / np.sqrt(g["n"])
        recs = []
        ref_base = g[(g["variant"] == "reference") & (g["role"] == "baseline")].set_index("scenario")
        for (v, sc), sub in g.groupby(["variant", "scenario"]):
            s = sub.set_index("role")
            if "candidate" not in s.index:
                continue
            bb = s.loc["baseline"] if "baseline" in s.index else ref_base.loc[sc]
            adv = bb["r"] - s.loc["candidate", "r"]
            se = math.sqrt(bb["se"] ** 2 + s.loc["candidate", "se"] ** 2)
            recs.append({"variant": v, "scenario": sc, "advantage": adv, "se": se, "advantage_rel": adv / bb["r"],
                         "t_ratio": adv / se if se > 0 else float("inf"),
                         "baseline_from": "same variant" if "baseline" in s.index else "reference"})
        adv = pd.DataFrame(recs)
        self.save(adv, "questions", "q3_advantage.parquet")
        ref = adv[adv["variant"] == "reference"].set_index("scenario")["advantage"]
        adv["relative_to_reference"] = adv.apply(lambda r: r["advantage"] / ref.get(r["scenario"], np.nan), axis=1)
        fpath = self.out / "failures.parquet"
        fails = {}
        if fpath.exists():
            f = pd.read_parquet(fpath)
            f = f[f["tag"].str.startswith("q3-")]
            fails = f.groupby("tag")["reason"].first().to_dict()
        res = {
            "advantage_by_variant": adv.groupby("variant")[["advantage_rel", "relative_to_reference", "t_ratio"]].mean().to_dict(orient="index"),
            "min_t_ratio_reference": float(adv[adv["variant"] == "reference"]["t_ratio"].min()),
            "failures": fails,
        }
        self.save(adv, "questions", "q3_advantage.parquet")
        _json(self.path("questions", "q3.json"), res)
        return res

    # ------------------------------------------------------------------ gate

    def gate(self) -> dict:
        ver = json.loads((self.out.parent / "verification" / "verification.json").read_text())
        cv = json.loads((self.out / "surrogate" / "cv_loop.json").read_text())
        comp = json.loads((self.out / "search" / "comparison.json").read_text())
        fin = json.loads((self.out / "final" / "final.json").read_text())
        reg = (self.root / "data" / "source_register.yaml").read_text()
        import yaml

        entries = yaml.safe_load(reg)["entries"]
        bn_ready = all(e["status"] in ("measured", "inferred") for e in entries if e["id"] in ("SR-001", "SR-009"))
        # compare methods on the confirmation batch: each method's best-of-four (chosen on the
        # verification batch) re-run on independent seeds. One campaign per method, so these
        # are statements about this campaign, not about the methods in general.
        mc = fin["method_comparison_confirmed"]
        p_adj = holm({k: v["p"] for k, v in mc.items()})
        resolved = {k: bool(mc[k]["diff"] > 0 and p_adj[k] < 0.05) for k in mc}
        checks = {
            "solver_checks_pass": all(r["passed"] for r in ver),
            "source_supported_bn_case_defined": bn_ready,
            "heldout_median_error_below_5pct": cv["gaussian_process"]["median_rel_error"] < TARGET_HELDOUT_ERROR,
            "finalist_top1_survives_direct_check": fin["top1_survives"],
            "loop_top1_not_resolved_from_verified_best": fin["top1_not_resolved_from_verified_best"],
            "bo_resolved_better_than_random_this_campaign": resolved["random_search_minus_bayesian_optimisation"],
            "bo_resolved_better_than_sweep_this_campaign": resolved["sweep_minus_bayesian_optimisation"],
            "numerical_error_below_claimed_benefit": fin["numerical_error_below_claim"],
        }
        required = ("solver_checks_pass", "source_supported_bn_case_defined", "heldout_median_error_below_5pct",
                    "finalist_top1_survives_direct_check", "bo_resolved_better_than_random_this_campaign",
                    "numerical_error_below_claimed_benefit")
        proceed = all(checks[k] for k in required)
        res = {
            "checks": checks,
            "proceed_to_phase2": proceed,
            "ml_earns_its_complexity": checks["bo_resolved_better_than_random_this_campaign"]
            and checks["bo_resolved_better_than_sweep_this_campaign"],
            "required": list(required),
            "blocking": [k for k in required if not checks[k]],
            "informational_unmet": [k for k, ok in checks.items() if not ok and k not in required],
            "method_picks": fin["method_picks"],
            "confirmed_pick_robust": {m: fin["confirmed"][l]["robust"] for m, l in fin["method_picks"].items()},
            "confirmed_pick_robust_se": {m: fin["confirmed"][l]["robust_se"] for m, l in fin["method_picks"].items()},
            "method_comparison_p_holm": p_adj,
            "solver_calls_to_within_1pct": {k: comp[k]["calls_to_within_1pct_of_overall_best"]
                                            for k in ("bayesian_optimisation", "random_search", "sweep")},
            "note": "single campaign per method; 'not resolved' is not equivalence",
        }
        _json(self.path("gate.json"), res)
        return res

    def manifest(self):
        write_manifest(self.out, self.config_dict())


def _time_generation(args) -> float:
    """CPU seconds to (re)run one geometry generation attempt."""
    from .geometry import generate, validate

    d, layer, gs, seed = args
    t0 = time.perf_counter()
    g = generate(d, layer, gs, seed)
    if g.ok:
        validate(g)
    return time.perf_counter() - t0


def _variance_reduction(df: pd.DataFrame, col: str, inputs: dict) -> dict:
    """Share of the spread of ``col`` removed by knowing each input exactly (grid average)."""
    total = float(df[col].var(ddof=0))
    out = {"total_sd": math.sqrt(total), "mean": float(df[col].mean()), "range": float(df[col].max() - df[col].min()),
           "fraction_positive": float((df[col] > 0).mean())}
    for name, c in inputs.items():
        if df[c].nunique() < 2:  # input not varied in the available data (e.g. infeasible levels)
            out[f"variance_reduction_{name}"] = float("nan")
            continue
        cond = float(df.groupby(c)[col].var(ddof=0).mean())
        out[f"variance_reduction_{name}"] = 1.0 - cond / total if total > 0 else 0.0
    return out
