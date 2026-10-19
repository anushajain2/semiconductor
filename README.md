# tcid: thermal composite inverse design

A research tool that recommends particle **size** and **orientation distributions** for a
ceramic-filled epoxy adhesive layer, then checks whether the thermal benefit survives
uncertain interfaces and geometric variability. The application is semiconductor packaging
and electronics thermal management. The tool designs a composite *microstructure*; it does
not invent a semiconductor or a polymer chemistry.

**Status: Phase 1 complete (2D verified solver and constrained design loop).** Phase 2
(3D oblate BN platelets with the anisotropic BN tensor, a source-supported BN material
case) has not been started. Phase 1 uses an isotropic *alumina-like development filler*
only to build and verify the method, so **none of its results is evidence about boron
nitride**. See [`docs/report.md`](docs/report.md) for the technical report and
[`docs/worked_example.md`](docs/worked_example.md) for how the output changes a design decision.

## What it does

| Component | Implementation |
|---|---|
| Geometry | Non-overlapping ellipses inside a finite layer with a resin-rich margin; five descriptors (loading, median diameter, log-normal width, aspect ratio, orientation order about a preferred direction). Impossible requests are rejected with a reason; achieved loading, minimum gap, seed and generator version are recorded. |
| Mesh | Gmsh conforming triangulation; polygon vertices on the exact ellipses; ≥ 3 elements across every narrow gap and ≥ 2 across the margin; automatic detection of single-element gap bridges. |
| Physics | Steady anisotropic conduction (scikit-fem P1 + SciPy). Particle–matrix interface resistance via duplicated interface DOFs and an interface-conductance term; perfect contact by tying DOFs; rotated particle tensors; isothermal faces, insulated sidewalls; outer contact resistances in series or as Robin conditions. |
| Verification | Homogeneous and layered analytic cases, energy balance, mesh refinement, rotation invariance, geometry validity, pure-phase limits, monotonicity, Wiener bounds (`tcid verify`). |
| Surrogate | Gaussian process over the design descriptors and uncertain physical inputs, predicting log resistance with seed-variance noise; extra-trees alternative; design-grouped cross-validation. |
| Search | Robust (worst-scenario) Bayesian optimisation with exploration points, compared against matched-budget random search and a 3-level sweep; shortlist re-simulated on fresh seeds and refined meshes. |
| Studies | Q1 design map vs interface resistance and orientation variability; Q2 contact-resistance decision boundary at three thicknesses; Q3 geometric perturbations; Q4 value of information; Phase 2 gate. |
| Traceability | Every solve stores its configuration hash, seed, achieved loading, particle statistics, properties, interface model, boundary conditions, mesh size, residual, heat balance, flux, resistance, timing, code version and status; cache keyed by the complete configuration; manifests with checksums. |

## Install

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e .[test]            # or: pip install -r requirements-lock.txt && pip install -e . --no-deps
```

The `gmsh` wheel loads system OpenGL/X libraries even when headless. On Debian/Ubuntu:
`apt-get install libglu1-mesa libxcursor1 libxinerama1 libxft2`.
Tested with Python 3.13 (exact versions in `requirements-lock.txt`).

## Use

```bash
tcid verify                                   # solver acceptance gates
tcid geometry --d50 10um --aspect 3 --order 0.8 --plot geom.png
tcid solve examples/case.yaml                 # one configured case over the scenario set
tcid all                                      # full Phase 1 pipeline (about 1 h on 4 cores)
tcid recommend --thickness "100 um" --loading 0.30
```

Pipeline stages can be run one at a time:
`tcid pilot | initial | search | final | questions | gate | figures | report`.
Each stage reads earlier outputs from `results/phase1/` and the solve cache in `cache/`.

## Layout

```
configs/            study, material case and scenario YAML
data/               source register (provenance of every physical range)
src/tcid/           geometry, mesh, solver, verification, surrogate, search, studies, figures, report, recommend, cli
tests/              pytest suite (analytic, conservation, rotation, geometry validity)
results/            verification, pilot, design loop, final verification, research questions (Parquet/JSON), manifest
docs/               report, worked example, pre-registered hypotheses, figures
examples/           a single-case configuration for `tcid solve`
```

## Limits

2D only; area fraction is not volume fraction; finite 200 µm-wide domain; all material
properties are *assumed* (see `data/source_register.yaml`). Stress scenarios are not
probabilities. Not modelled: processing, adhesion, cure, dielectric breakdown, stress,
ageing. This is not a qualified adhesive and not validated semiconductor reliability.
