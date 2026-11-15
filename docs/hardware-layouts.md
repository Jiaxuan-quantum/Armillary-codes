# Hardware layouts of the three representative codes

The final layouts preserve the code matrices, all physical pairings, and the forward/reverse syndrome-extraction circuit. Geometry changes neither the candidate definitions nor the published circuit-level noise model. The independent module in [`hardware_layouts/`](../hardware_layouts/) supplies the continuous routes and their allocation to routing faces and tiers.

## Results

All lengths and clearances below use the minimum center-to-center distance `a0` among all `2n` physical qubit positions. Each layout has `a0 = 1`.

| Code | Physical connections | Mean length / a0 | Maximum length / a0 | Wire separation lower bound / a0 | Device clearance lower bound / a0 |
|---|---:|---:|---:|---:|---:|
| `[[3636,1822,12]]` | 24,038 | 9.3237 | 18.8002 | 0.100708 | 0.097228 |
| `[[4932,2470,14]]` | 32,195 | 10.3428 | 19.0539 | 0.100539 | 0.099916 |
| `[[8688,4348,16]]` | 66,065 | 18.8947 | 48.4811 | 0.100097 | 0.065302 |

The lower bounds are rounded down. The mean in this table weights all physical connections equally. The tier-averaged Length, tiers, TSVs, Bonds, and combined score are displayed in the [main README](../README.md#hardware-layouts). Machine-readable reference values are in [`reference_results.json`](../hardware_layouts/reference_results.json).

## Reproduce the layouts

With Python 3.11 or newer and the repository's base dependencies installed, run:

```bash
python hardware_layouts/reproduce.py --out results/layouts
python hardware_layouts/reproduce.py --code d12 --out results/layout_d12
python hardware_layouts/reproduce.py --code d14 --out results/layout_d14
python hardware_layouts/reproduce.py --code d16 --out results/layout_d16
```

The module can also be copied to a separate directory and run with NumPy and SciPy installed. Use a new output directory for each run. The selected layout inputs are immutable, and regenerated results are written to the requested output directory. Output includes `SUMMARY.json`, and each code's `METRICS.json`, `VERIFY.json`, and reconstructed `routes.json`.

Verification checks file checksums, reconstructs the code matrices and their identity, replays both circuit directions, checks the complete physical connection set, regenerates continuous routes, verifies geometric separation, and recomputes resource metrics. It also enumerates 840 possible stackings for each eight-face layout and 10,395 for the eleven-face layout. Saved logical witnesses are checked. Archived distance lower bounds are inherited by matrix identity, rather than independently reproved.

## Geometry and routing data

Physical positions use the manuscript coordinates `(r,c,i,j)`. The module with label `c` contains `(0,c,i,j)` and `(1,c-1,i,j)`. For `M = 12`, the inner-to-outer module order is `(0,11,1,10,2,9,3,8,4,7,5,6)`. Channels increase outwards in `j`, and the order of the two rails reverses across the fold.

The data store one representative of each connection under station translation. Its `p` rotated copies generate all physical connections. The fixed-edge set includes SWAP-only connections associated with unretained check groups. Repeated uses of the same pair count as one physical connection.

The continuous path is specified by consecutive vertices in unfolded coordinates `(u,R)`, with `u = R0*theta`. Each straight segment in this coordinate system is mapped continuously to the annulus by

```text
(u, R) -> (R cos(u/R0), R sin(u/R0)).
```

Consequently, a segment with fixed `R` represents a circular arc, and a segment changing both coordinates is generally curved in the physical plane. The physical lengths include endpoint leads and integrate the mapped curve. Treating the exported vertices as straight Cartesian chords would give different lengths and clearance tests.

Geometry data give local radial spacings and angular offsets. Hardware data bind each physical connection to a face. Radial and circumferential paths can share a face, and every path remains on that face throughout its in-plane route. All station copies use the same face. Each endpoint has at most four connections on any one face.

One face occupies the qubit tier. Each higher tier can host two faces, with an empty face where needed. For a fixed pairing, the more populated face is placed below the other face, and more populated tiers are placed closer to the qubit tier. The verifier enumerates the choice of bottom face and pairing to confirm the minimum HAL score for these fixed face assignments. This is not a globally optimal solution over geometries and face assignments.

## Continuous clearance

Each physical device has a fixed region defined in unfolded coordinates by `|u-u_c| <= d_u` and `|R-R_c| <= d_r`, then mapped to the annulus. The half-widths `(d_r,d_u)/a0` are `(0.152,0.405)`, `(0.152,0.3675)`, and `(0.35,1.005)` for distances 12, 14, and 16. Since `d_u` is measured at the reference radius, its physical arc extent varies with radius.

Two paths on the same face may meet only inside their common endpoint regions. Outside those regions their separation must exceed `0.1*a0`. Every complete path must stay more than `0.05*a0` from every device region other than its own two endpoints. Only shared endpoint regions are exempted from the path-pair check. Nonendpoint device checks include the full path and its leads.

The audits combine exact line-segment calculations in unfolded coordinates with analytic bounds for the polar map and periodic translations. They check continuous curves, not just sampled vertices. Calculations use double precision and are not formal interval-arithmetic certificates. The reported bounds describe centerlines and fixed device regions, without manufacturing tolerances, shielding, vertical-interconnect footprints, or control/readout wiring.

## HAL metrics

`Length` is the equally weighted mean of the mean connection lengths in each occupied tier, including the qubit tier. Each physical path is counted once. Because tiers are equally weighted, moving unchanged routes between tiers can change Length without changing the all-connection mean.

`TSVs` is the mean TSV count over connections routed above the qubit tier. `Bonds` is the maximum tier mean of HAL's face-switch segment count. It is not the number of physical bump bonds. Avoiding intermediate face switches does not eliminate the vertical access needed to reach upper tiers or their upper surfaces. The output also separates total TSV count, physical face switches, and full-stack bonds under the stated stack-access convention.

```text
C_hw = 1 + (T - 1)/16 + max(0, Length - 1)/36 + TSVs/12 + Bonds/16
```

The tier aggregation follows HAL, while the length estimator uses complete continuous arc length normalized by `a0`. Published HAL routing counts axial and diagonal grid steps equally. These length estimators should not be described as identical.

## Provenance and scope

The release records the original frozen input hashes and any format adaptations in `hardware_layouts/PROVENANCE.json`. Its manifest protects the distributed reproduction inputs and implementation. The code identities agree with the 122-candidate catalog. Historical search logs, obsolete layouts, solver installation details, and private filesystem paths are not needed to reproduce the selected layouts and are not distributed.

These are geometric layouts with numerical continuous-separation checks. Reproduction does not run a new layout search, an exhaustive distance proof, a memory experiment, or an electromagnetic simulation. The original QEC results use the stated uniform circuit noise model and do not predict length-dependent hardware error rates.

## Reference

The hardware architecture and aggregate score follow M. Mathews et al., *Placing and routing quantum LDPC codes in multilayer superconducting hardware*, npj Quantum Information **12**, 114 (2026), [doi:10.1038/s41534-026-01243-w](https://doi.org/10.1038/s41534-026-01243-w), and the [HAL implementation](https://github.com/EQuS/HAL). The Armillary layout verifier independently evaluates its continuous routes and does not invoke the HAL place-and-route search.
