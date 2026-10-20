| Check | Proposed acceptance criterion | Worst value | Pass |
|---|---|---|---|
| homogeneous_layer | R = t/k within 0.5% | 1.76e-15 | yes |
| layered_interfaces | summed layer + interface resistances within 1% | 1.01e-13 | yes |
| energy_balance | |Q_in - Q_out| / Q < 0.1% | 8.64e-12 | yes |
| mesh_refinement | resistance change < 2% on the next refinement | 2.50e-04 | yes |
| orientation_implementation | rigid rotation exact (<1e-8); 90-degree re-meshed rotation within 1%; un-rotated tensor detected | 1.26e-05 | yes |
| geometry_validity | no overlap or margin violation in accepted geometries; achieved loading within tolerance and recorded | 0.00e+00 | yes |
| limits_and_bounds | pure-phase exact; monotone in R_interface; R->inf matches holes within 1%; inside Wiener bounds | 1.34e-14 | yes |
