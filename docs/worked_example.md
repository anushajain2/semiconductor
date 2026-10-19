# Worked example: how the output changes a material-design decision

> Development case only: a 2D, alumina-like isotropic filler stand-in. The *form* of the
> decision carries over to BN; the numbers do not.

## The situation

A packaging group bonds a power module to a copper spreader with a ceramic-filled epoxy.
The bond line is 100 µm and the formulation is fixed at 30% filler (area fraction in this 2D
model). The filler supplier offers several particle grades, and the process team can add a
field-assisted alignment step. That step costs money and may only partly align the
particles. The group wants to know three things:

1. Is alignment worth paying for?
2. Which particle grade should they buy?
3. Is any of this worth doing, given that the two outer contacts also add resistance?

## Ask the tool

```bash
tcid recommend --material alumina_epoxy_dev_2d --thickness "100 um" --loading 0.30
```

```
candidate                                 d50 um  sigma    AR     S       R_robust    R_nom  vs baseline
nominal_opt *                               6.72   0.42  4.00  0.87    180.8mm2K/W    167.5      +22.4%
robust_opt ~                               21.09   0.05  4.00  0.84    186.9mm2K/W    179.4      +19.7%
shortlist_2 ~                               9.61   0.40  4.00  0.95    183.6mm2K/W    172.7      +21.2%
shortlist_3                                16.65   0.17  4.00  0.95    185.5mm2K/W    177.1      +20.3%
baseline (random orientation)               6.72   0.42  4.00  0.00    232.8mm2K/W    215.7       +0.0%

* recommended;  ~ statistically indistinguishable from the recommended design (< 2 SE): choose on other grounds.
```

The tool also writes `results/recommendation/*.yaml`. These are configuration files that
re-run any candidate with `tcid solve`. It also writes a CSV of per-scenario resistances and
the statement of the modelled domain.

## What the output changes

**Decision 1, alignment: yes, and it is the decision that matters.** Every aligned,
elongated candidate (aspect ratio 4, orientation order 0.84–0.95) cuts the worst-case layer
resistance by 19–23% relative to the same particles randomly oriented. That holds in every
interface and filler scenario and survives fresh realisations (t = 5.6–16). The perturbation study
(Q3) puts a price on imperfect processing: if the achieved orientation order is 0.2 lower than
targeted, about three quarters of the benefit remains. The design map (Q1) shows that random
orientation never wins anywhere in the scenario space. Without the tool, a team might have
debated grades first. The tool says the alignment step is where the thermal value is.

**Decision 2, grade: do not pay a premium for a particular size yet.** The recommended
6.7 µm grade and the 9.6 µm, 16.6 µm and 21 µm alternatives lie within about 3% of each other.
Two of them are statistically indistinguishable from the recommendation. The design map
explains why the choice is unresolved. Smaller particles win if the particle–matrix interface
is good (below about 6e-7 m²K/W). Larger particles win if it is poor (above about 1.5e-6).
The value-of-information analysis (Q4) puts 92% of the size decision on interface resistance.
The actionable output is therefore a measurement, not a grade: characterise the effective
interface resistance of the candidate grades before choosing between them. Until then, choose
the grade on cost, viscosity or supply.

**Decision 3, is it worth it at all?** With identical outer contacts the ranking cannot
change, only the size of the benefit. For this 100 µm layer the improvement in assembly
resistance stays above 10% until the two contacts together exceed about 330 mm²K/W. That is
far above the assumed range for well-made contacts (up to about 60 mm²K/W). For a 50 µm bond
line the threshold falls to about 99 mm²K/W. A thin-bond-line product with poor contacts
should fix its contacts first.

## Outside the verified domain

```bash
tcid recommend --thickness "60 um" --loading 0.25 --run
```

The tool answers `UNSUPPORTED REQUEST` and lists the reasons: the thickness and loading are
outside the verified set. With `--run` it also solves the verified candidates at the new
conditions, clearly labelled *not re-optimised*. Here the recommended distribution still
beats its random-orientation baseline (134 vs 159 mm²K/W worst case, 3 seeds). That is a
direct solver result, not a validated recommendation for that layer.

## What this example does not say

It says nothing about BN, about 3D connectivity, about how alignment is achieved or what it
costs, or about adhesion, cure, dielectric strength or reliability. Phase 2 must repeat the
decision with 3D oriented BN platelets and a source-supported material case before any
number above can inform a real formulation.
