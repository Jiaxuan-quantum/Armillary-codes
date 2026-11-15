# Third-party acknowledgments

This project uses the following software as dependencies. Their code and binary
distributions are not vendored here, and their original licenses continue to apply.

| Software | Use | License / source |
|---|---|---|
| NumPy | Numerical arrays | BSD-3-Clause; https://numpy.org |
| SciPy / HiGHS | Sparse matrices, statistics, and open-source integer optimization | SciPy: BSD-3-Clause; HiGHS: MIT; https://scipy.org and https://highs.dev |
| Stim | Stabilizer circuits and detector error models | Apache-2.0; https://github.com/quantumlib/Stim |
| Relay-BP 0.2.2 | First-stage belief propagation | Apache-2.0; https://github.com/trmue/relay |
| COPT / coptpy (optional) | Integer optimization | Solver license supplied separately by Cardinal Operations; https://www.shanshu.ai/copt |
| Matplotlib (optional) | Benchmark plotting | Matplotlib license; https://matplotlib.org |
| pytest (development) | Tests | MIT; https://pytest.org |

Relay-BP is authored by Thomas Alexander and Tristan Müller. Its underlying method
is described by Tristan Müller and collaborators in *Improved Belief Propagation
Is Sufficient for Real-Time Decoding of Quantum Memory*, arXiv:2506.01779 (2025).
Stim is described by Craig Gidney in *Stim: a fast stabilizer circuit simulator*,
Quantum **5**, 497 (2021), https://doi.org/10.22331/q-2021-07-06-497.

The Armillary construction, circuit adapter, experiment interface, and curated
candidate data in this repository are project materials covered by its LICENSE.

The standalone hardware-layout module uses the architecture and score conventions of M. Mathews et al., *Placing and routing quantum LDPC codes in multilayer superconducting hardware*, npj Quantum Information **12**, 114 (2026), https://doi.org/10.1038/s41534-026-01243-w. The reference HAL implementation is https://github.com/EQuS/HAL. HAL is not a runtime dependency of the layout verifier, and its router is not vendored here. The continuous-path and layout-verification code is original project code under this repository's MIT license.
