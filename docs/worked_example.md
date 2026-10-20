# Worked example: how the output changes a material-design decision

> Development case only: a 2D, alumina-like isotropic filler stand-in. The *form* of the
> decision carries over to boron nitride (BN); the numbers do not.

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
candidate                           d50 um  sigma    AR     S   verif.       confirmed vs ref. base
sweep_3 *                            24.00   0.05  4.00  0.95    175.4    178.5 ±  3.0      +20.2%
sweep_4 ~                             9.80   0.05  4.00  0.95    179.8    181.2 ±  0.9      +18.9%
nominal_opt ~                         6.72   0.42  4.00  0.87    180.8    182.9 ±  0.8      +18.2%
sweep_2 ~                             9.80   0.50  4.00  0.95    181.7    179.9 ±  1.4      +19.5%
sweep_best ~                          9.80   0.28  4.00  0.95    183.1    179.3 ±  0.9      +19.8%
shortlist_2 ~                         9.61   0.40  4.00  0.95    183.6    180.8 ±  2.1      +19.1%
shortlist_3 ~                        16.65   0.17  4.00  0.95    185.5    179.2 ±  1.5      +19.9%
robust_opt ~                         21.09   0.05  4.00  0.84    186.9    185.0 ±  2.2      +17.3%
random_search_2                      11.20   0.07  3.40  0.91    193.9    191.4 ±  1.0      +14.4%
baseline (random orientation)        24.00   0.05  4.00  0.00    237.2    223.6 ±  3.3       +0.0%

* recommended;  ~ difference from the recommended design not resolved
```

How to read the table:

- Values are worst-scenario layer resistance in mm²K/W.
- "verif." is the 6-seed batch the candidates were chosen on.
- "confirmed" is a separate 12-seed batch that was not used to choose anything. Quote these values.
- The last column compares each candidate with the recommended design's own random-orientation baseline.

The tool also writes `results/recommendation/*.yaml` (re-run any candidate with `tcid solve`),
a CSV of per-scenario resistances and a statement of the modelled domain.

## What the output changes

**Decision 1, alignment: yes. This is the decision the simulation supports firmly.**

- Each aligned, aspect-ratio-4 candidate has 19–27% lower worst-case resistance than the same particles randomly oriented.
- On the independent confirmation batch, the recommended design's advantage over its baseline is 45 mm²K/W (95% CI 36–54), about 20%.
- In the perturbation study, lowering the achieved orientation order by 0.2 did not reduce the advantage of this design.
- Random orientation never wins a cell of the design map.

**Decision 2, grade: the simulation cannot choose one. That is itself the useful answer.**

- The top row is only the mechanical pick. Every candidate marked `~` is statistically unresolved from it, and the confirmed values of the aspect-ratio-4 designs span only 178.5–185.0 mm²K/W.
- With no thermal basis for choosing, the tool says to choose on other grounds. Two of them come from the simulation itself:
  * **Thickness flexibility.** The 24 µm pick cannot be used in a 50 µm or 60 µm layer: its particles are longer than the size rules allow, as the run below shows. The ~10 µm designs fit thinner layers.
  * **Variability.** The 24 µm pick varies most from one realization to the next (± 3.0 vs about ± 0.9). Its advantage also drops to 0.44× if the particle size runs 20% larger than specified.
- A team that might use the adhesive in several bond-line thicknesses would reasonably buy a ~10 µm grade, such as `sweep_best`.

**When does particle size start to matter thermally?**

- Only if the particle–epoxy interface is poor. Below roughly 7–9e-7 m²K/W of interface resistance, smaller particles do slightly better; above it, larger particles do.
- That crossover is uncertain: its 95% intervals run from about 2e-7 to 2e-6, straddling the top of the supported range.
- So if the interface resistance of the candidate grades might be near or above 1e-6 m²K/W, characterising it is what would make size matter. Below that, the size choice has little thermal consequence.

**Decision 3, is it worth it at all?**

- With identical outer contacts the ranking cannot change, only the size of the benefit.
- For this 100 µm layer the assembly benefit is 20% [17–24%] with ideal contacts. It stays above 10% until the two contacts together reach about 220 mm²K/W (95% interval 144–299).
- Even the bottom of that interval is more than twice the assumed range for well-made contacts (up to about 60 mm²K/W).
- Thicker layers tolerate even more contact resistance. Thin bond lines could not be evaluated: the recommended distribution does not fit in 50 µm, and the generator failed for part of the medium design's realizations there.

## Outside the verified domain

```bash
tcid recommend --thickness "60 um" --loading 0.25 --run
```

The tool answers `UNSUPPORTED REQUEST`, because the thickness and loading are outside the verified set. With `--run` it also tries to solve the verified recommendation at the new conditions, clearly labelled *not re-optimised*. Here it reports that the recommended distribution cannot be placed at all at 60 µm under the size rules. That is a concrete reason not to transfer this recommendation to thinner bond lines without a new study.

## What this example does not say

It says nothing about:

- boron nitride;
- 3D connectivity;
- how alignment is achieved, or what it costs;
- adhesion, cure, dielectric strength or reliability.

Phase 2 must repeat the decision with 3D oriented BN platelets and a source-supported material case before any number above can inform a real formulation.
