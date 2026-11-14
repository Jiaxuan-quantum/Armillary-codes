# Two-stage decoding

The detector error model supplies a binary detector matrix `H` and independent
fault probabilities `p`. A decoder receives only these quantities and the
observed detector syndrome `s`. Logical observables and sampled logical outcomes
are used separately to score a returned correction.

## Relay-BP

The first stage uses the float64 decoder from
[`relay-bp` 0.2.2](https://github.com/trmue/relay). It runs a short Relay pass and
checks the returned correction explicitly against `H e = s (mod 2)`. A long
pass is run only if that check fails. If neither pass gives a valid correction,
the sample is handed to MLE decoding.

| Parameter | Short pass | Long pass |
|---|---:|---:|
| `pre_iter` | 200 | 500 |
| `num_sets` | 20 | 200 |
| `set_max_iter` | 100 | 200 |
| `gamma0` | 0.125 | 0.35 |
| `stop_nconv` | 2 | 2 |
| `seed` | 0 | 0 |

The same strategy is used for both memory bases of the three highlighted codes.
All other Relay parameters retain the defaults of version 0.2.2. Each pass is
constructed anew for each sample. Zero syndrome is assigned the zero correction.
A valid syndrome does not certify logical correctness: a converged but wrong
correction is a first-stage decoding failure and is not sent to the second stage.

Relay-BP was introduced by T. Müller et al., *Improved Belief Propagation Is
Sufficient for Real-Time Decoding of Quantum Memory*,
[arXiv:2506.01779](https://arxiv.org/abs/2506.01779). The upstream implementation
is distributed under Apache-2.0; it is an installed dependency, not vendored here.

## MLE

For an independent detector error model, the most likely error pattern minimizes

```math
\sum_j \log\!\left(\frac{1-p_j}{p_j}\right)e_j,
\qquad H e - 2z = s,\quad e_j\in\{0,1\},\quad z_i\in\mathbb Z.
```

This repository provides a portable SciPy/HiGHS implementation and an optional
COPT backend. A suitable COPT license is required for research-size models.
Install the `copt` extra and choose `backend="copt"` to use it. Every returned
correction is checked against the binary syndrome before use. Columns without
detector support are retained. The decoder solves a single error-pattern MLE
problem; it does not sum probabilities over degenerate logical classes.

```python
import numpy as np
from armillary.decoders import MLEDecoder, RelayDecoder

H = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.uint8)
probabilities = np.array([0.1, 0.1, 0.1])
syndrome = np.array([1, 1], dtype=np.uint8)

answer = RelayDecoder(H, probabilities).decode(syndrome)
if not answer["accepted"]:
    answer = MLEDecoder(H, probabilities, backend="scipy").decode(syndrome)
```

`accepted` means that a syndrome-consistent correction was returned. The MLE
diagnostics also identify the solver's termination status. Logical success must
be evaluated from the correction's predicted observables. Execution errors and
absent corrections are not counted as successful decodes.

This is a reference workflow for the published construction and decoding method.
The included aggregate measurements describe the research experiments; running
the reference workflow generates new samples rather than replaying those records.
