# Reproducible hardware layouts

This standalone module contains the final continuous wiring of the three representative Armillary codes. It is distinct from the 122-entry code catalog, whose discrete routing rules determine circuit pairings but do not specify wire geometry.

| Code | Tiers | Faces | HAL Length | TSVs | Bonds | C_hw |
|---|---:|---:|---:|---:|---:|---:|
| `[[3636,1822,12]]` | 5 | 8 | 7.6116 | 3.4200 | 0.9063 | 1.7753 |
| `[[4932,2470,14]]` | 5 | 8 | 8.5930 | 3.4300 | 1.0096 | 1.8098 |
| `[[8688,4348,16]]` | 6 | 11 | 15.7464 | 4.4255 | 1.2333 | 2.1680 |

Length uses complete continuous arc length, normalized by the minimum physical qubit pitch `a0`, with equal-weight tier averaging. It is distinct from the all-connection mean and from HAL's published equal-grid-step length estimator. Bonds is HAL's segment metric, not a physical bump-bond count.

## Run

From the repository root:

```bash
python -m pip install numpy scipy
python hardware_layouts/reproduce.py --out results/layouts
# One code only
python hardware_layouts/reproduce.py --code d16 --out results/layout_d16
```

The `hardware_layouts` directory can be copied elsewhere. With Python 3.11 or newer, NumPy, and SciPy installed, run `python /path/to/hardware_layouts/reproduce.py --out /path/to/results` from any working directory. No solver, network, Relay-BP, or original research tree is needed. Do not use Python's `-O` option, which disables assertions.

The command checks the manifest, reconstructs matrices and both circuit directions, rebuilds all routes including SWAP-only connections, checks continuous separation, and recomputes the resource metrics and the best stacking for the fixed face assignment. It writes `SUMMARY.json` and per-code `METRICS.json`, `VERIFY.json`, and `routes.json` to a new output directory.

Wire separation outside shared endpoint regions must exceed `0.1*a0`. Full paths must remain more than `0.05*a0` from nonendpoint device regions. The calculations certify numerical analytic lower bounds on continuous curves, rather than sampling vertices. They use double precision and are not formal interval-arithmetic proofs.

`reference_results.json` contains the frozen reference values. `PROVENANCE.json` records source hashes and format adaptations, while the manifest checks the distributed inputs and implementation. The data retain code identities, device geometry, physical connection representatives under station translation, continuous paths, face assignments, and tier stacking.

These checks reproduce the selected layouts, not the optimization search that found them. Fixed-face stackings are enumerated exhaustively, while the overall geometry and face assignments are not proven globally optimal. Saved logical witnesses are replayed, but distance lower bounds are inherited by code identity. This module does not rerun QEC simulations or model electromagnetic behavior, control/readout routing, via footprints, or manufacturing tolerances.

See the repository's [hardware-layout documentation](../docs/hardware-layouts.md) for coordinate and metric definitions. The original project [MIT license](../LICENSE) applies to this module.
