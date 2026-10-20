# Phase gate v2 and corrected-surrogate campaign plan (written before the runs they govern)

Written on 2026-10-20 after an external review, **before** the corrected-surrogate search
campaign below was run and before any Phase 2 work. This file is not edited after those runs.

## Why a new gate

The Phase 1 gate (results/phase1/gate.json) required the loop's top-ranked finalist to keep
first place on fresh seeds. With several candidates within noise of each other, that criterion
mostly tests luck. **Its failure stands and is not re-decided**: Phase 1 did not pass its gate.
The criteria below replace it for future evaluations (Phase 2 onward); they are not applied
to re-judge Phase 1.

## Gate v2 (applies from Phase 2)

A recommendation passes only if every condition holds, each evaluated on a confirmation batch
of at least 12 independent realizations that was not used to choose anything:

| ID | Condition |
|---|---|
| G1 | All solver verification checks pass at the code version that produced the results (clean, tagged commit). |
| G2 | The material case is source-supported: the register entries it uses are `measured` or `inferred`, not `assumed` or `pending`. |
| G3 | Retained advantage: the lower bound of the 95% CI of (matched random-orientation baseline − recommended) / baseline, worst stress scenario, is at least 10%. |
| G4 | Acceptable standing among alternatives: no confirmed alternative is resolved better than the recommended design (Holm-corrected Welch, p < 0.05), and the upper 95% CI bound of (recommended − best confirmed alternative) is at most 3% of the recommended value. |
| G5 | Numerical error (next mesh refinement) is below 20% of the smallest benefit claimed. |
| G6 | If a surrogate guided the search: design-grouped held-out median error < 5% **and** 95% interval coverage between 85% and 99%. |

Re-running a campaign with new seeds or settings until a preferred design passes is not
allowed; every run that counts toward a gate is recorded, including failed ones.

## Corrected-surrogate campaign (methods follow-up, Phase 1 data)

* Defect: per-point observation noise was supplied to a GP fitted with `normalize_y=True`
  without dividing by the target variance (≈39× too small); the fitted white-noise term
  collapsed to its lower bound.
* Plan: fix the scaling; re-run design-grouped cross-validation (the uncorrected results are
  kept as `*_v1_uncorrected`); run **one** new BO campaign (`bo_v2`) from the same 40-design
  initial set, same budget (10 iterations × 4 designs × 3 seeds), same seeds and acquisition
  settings as the original campaign.
* Evaluation: its top four designs by loop robust value are verified on the existing
  verification seeds (1000–1005); its best-of-four is confirmed on the existing confirmation
  seeds (6000–6011) and compared with the existing confirmed picks.
* The Phase 1 recommendation and final selection are **frozen**: the new campaign is reported
  as a methods comparison and does not re-open the choice of recommended design.
