# Armillary codes

Reproducible definitions, circuit-level memory experiments, and reference decoders for **rate-1/2 quantum LDPC codes for multi-tier superconducting hardware**.

This repository accompanies the Armillary-code study. It contains a compact catalog of **122 candidates**, a physical-site CXSWAP circuit builder, short-to-long Relay-BP decoding, a most-likely-error (MLE) decoder, and independently reproducible hardware layouts for the three representative codes. The three representative codes are:

| Identifier | Code parameters `[[n, k, d]]` | Catalog identifier |
|---|---|---|
| `8749` | `[[3636, 1822, 12]]` | `P-2b28b13c8749` |
| `9053` | `[[4932, 2470, 14]]` | `P-54f8f55d9053` |
| `p181` | `[[8688, 4348, 16]]` | `P-0bc3b705dfd0` |

## Install

Use Python 3.11 or newer in a virtual environment:

```bash
git clone https://github.com/Jiaxuan-quantum/Armillary-codes.git
cd Armillary-codes
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[relay,test]"
```

On Windows, activate with `.venv\Scripts\Activate.ps1` in PowerShell. Relay-BP is a native dependency, pinned to version **0.2.2**. Wheel availability depends on Python and platform; installing from source requires a Rust toolchain. Linux x86-64 with Python 3.14 and glibc 2.34 or newer has a published wheel. The base installation, `pip install -e .`, supports candidate reconstruction, circuit generation, statistics, and SciPy/HiGHS MLE without Relay-BP.

For the optional **COPT** MLE backend:

```bash
python -m pip install -e ".[copt]"
```

Configure your COPT license through the solver's usual installation procedure. No solver license is distributed with this project. SciPy/HiGHS requires no commercial license.

## Quick start

### 1. Inspect and reconstruct the candidates

```bash
armillary catalog --distance 12
armillary verify 8749
armillary verify
```

The last command reconstructs all 122 candidates, checks the CSS matrices and their canonical SHA-256 values, and verifies the forward/reverse routing. Catalog coverage is 18 distance-12, 96 distance-14, and 8 distance-16 candidates. Reported distances are archived construction results; this command does not rerun distance searches. See [candidate definitions](docs/candidates.md).

### 2. Build a circuit and sample its first-stage decoder

```bash
armillary model 8749 --basis x --p 0.001 --rounds 4 --out work/8749_x
armillary sample work/8749_x --shots 32 --seed 1001 --out results/8749_x
```

The first command saves a Stim circuit and its detector model. Model construction can take a few minutes and uses more memory for larger candidates. The second command samples the **circuit**, runs Relay-BP, and reports correct outputs, incorrect outputs, and abstentions. An abstention is a handoff to the MLE stage, not a logical failure by itself.

For a matching Z-memory experiment, build with `--basis z`, use the same `--shots`, a different seed (for example `1002`), and separate model/result directories. Keep a new seed and directory for each independent run. Sampling streams also depend on Stim version, CPU architecture, and batch size.

### 3. Decode the saved handoffs

```bash
# Open-source reference backend
armillary decode work/8749_x results/8749_x --backend scipy --out results/8749_x_mle

# Or use COPT on the same saved input
armillary decode work/8749_x results/8749_x --backend copt --out results/8749_x_copt
```

This stage minimizes the negative log-likelihood of an independent detector-error pattern subject to the observed syndrome. Large instances can be computationally demanding. Both backends independently verify the returned binary correction against the syndrome. Simulated logical outcomes are used only for scoring; neither decoder receives them. See [decoding methods](docs/decoding.md).

### 4. Inspect the outputs

| File | Contents |
|---|---|
| Model directory | `model.stim`, detector matrix, fault priors, packed logical map, and checksummed metadata |
| `results/8749_x/summary.json` | Completed first-stage counts, requested remainder, seed, and model identity |
| `results/8749_x/handoffs.npz` | Syndromes and sample identifiers for Relay abstentions |
| `results/8749_x/scoring.npz` | Separate simulated logical outcomes for those handoffs |
| MLE output directory | Per-sample corrections, outcomes, and a combined summary when all requested samples are resolved |

Output directories must be empty. Interruptions and execution errors are reported separately from decoded outcomes. The reference commands are single-process examples; they do not launch cluster jobs.

## Methods and benchmark data

The circuit model uses four noisy QEC rounds with ideal boundaries by default. Single-qubit depolarization acts on ancillary preparation/readout and once on each data qubit per measurement/preparation interval. Two-qubit depolarization acts after each complete CXSWAP or SWAP-only operation. Forward and reverse routing alternate. See [circuit conventions](docs/circuits.md).

The front decoder first runs a short Relay-BP configuration. It accepts the first syndrome-consistent correction, otherwise runs a fresh longer configuration with the original priors. If both attempts fail the syndrome check, MLE receives the sample. The same two Relay configurations are used for X and Z memory experiments of all three representative codes.

The released aggregate data at physical error probability `p = 0.001` are:

| Code | Completed X/Z samples | Final block failures | Logical error rate per logical qubit per round |
|---|---:|---:|---:|
| `8749` | 6,723,292 | 22 | 4.49 × 10⁻¹⁰ |
| `9053` | 20,773,918 | 9 | 4.38 × 10⁻¹¹ |
| `p181` | 35,484,090 | 1 | 1.62 × 10⁻¹² |

We report `P_L = 1 - (1 - F/N)^(1/(k r))`, with `r = 4`. The distance-16 estimate is based on one failure and does not establish a probability below 10⁻¹². See [benchmark statistics](docs/results.md) for per-basis counts and uncertainty.

```bash
armillary benchmark
python examples/plot_benchmark.py --out results/benchmark.pdf
```

Install `pip install -e ".[plot]"` for the plotting example. The benchmark command recalculates the reported estimates from packaged counts; it does not rerun the sampling campaign. This release reproduces the code definitions and reference workflow, rather than every numerical decision in the original computation.

## Hardware layouts

The [hardware-layout reproduction module](hardware_layouts/) separately provides the final continuous wiring for the three representative codes. The 122 candidate definitions specify code construction and circuit pairings. Complete geometric layouts are released for these three codes only.

| Code `[[n, k, d]]` | Tiers | Routing faces | Length | TSVs | Bonds | C_hw |
|---|---:|---:|---:|---:|---:|---:|
| `[[3636, 1822, 12]]` | 5 | 8 | 7.6116 | 3.4200 | 0.9063 | **1.7753** |
| `[[4932, 2470, 14]]` | 5 | 8 | 8.5930 | 3.4300 | 1.0096 | **1.8098** |
| `[[8688, 4348, 16]]` | 6 | 11 | 15.7464 | 4.4255 | 1.2333 | **2.1680** |

Length follows HAL's tier averaging, using **complete continuous arc length divided by the minimum qubit pitch**. It is distinct from HAL's published equal-grid-step estimator and from the mean over all connections. TSVs and Bonds retain HAL's aggregation rules. Bonds is a source-defined segment metric, not a physical bump-bond count.

From a repository checkout, with the base dependencies installed:

```bash
python hardware_layouts/reproduce.py --out results/layouts
# Or reproduce one layout
python hardware_layouts/reproduce.py --code d16 --out results/layout_d16
```

This reconstructs the complete paths, checks circuit-edge coverage including SWAP-only connections, verifies continuous wire and device clearance, and recalculates the HAL metrics and optimal stacking for the fixed face assignment. The requirements are wire separation greater than `0.1` and nonendpoint-device clearance greater than `0.05`, in units of the minimum qubit pitch. No solver license, Relay-BP installation, network access, or original research directory is needed.

The module is a standalone part of this source repository and is not included in the Python wheel. It reproduces the selected layouts, without rerunning layout searches, distance certification, or QEC sampling. See [layout data and methods](docs/hardware-layouts.md) for coordinates, file formats, numerical checks, and limitations.

## Project layout

```text
src/armillary/
  data/             Candidate definitions, catalog, aggregate benchmark counts
  algebra.py        Binary linear algebra and logical bases
  candidates.py     Literal routing replay and CSS reconstruction
  circuits.py       Physical-site Stim memory circuits
  models.py         Detector matrices, priors, and logical maps
  decoders.py       Relay-BP and MLE reference implementations
  experiments.py    Sampling, handoffs, and independent scoring
  statistics.py     Failure-rate estimates
hardware_layouts/   Three final continuous layouts and standalone reproduction
docs/               Construction, circuit, decoding, statistics, and layout conventions
examples/           Small runnable examples
tests/              Construction, mapping, decoding, and accounting checks
```

Run `python -m pytest -q` for the test suite. See [validation](docs/validation.md) for the checked scope.

## Attribution and license

Original project code is provided under the [MIT license](LICENSE). Dependencies retain their own licenses. In particular, Relay-BP is the implementation by Thomas Alexander and Tristan Müller, available at [trmue/relay](https://github.com/trmue/relay), and Stim is Craig Gidney's stabilizer-circuit simulator. See [third-party acknowledgments](THIRD_PARTY.md).

When using this repository, cite the accompanying Armillary-code study and the simulation/decoding software relevant to your experiment. Repository citation metadata is in [CITATION.cff](CITATION.cff).
