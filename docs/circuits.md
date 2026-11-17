# Circuit simulation and detector models

The memory experiment is generated in fixed physical-site coordinates. Data and
auxiliary roles are tracked through the actual SWAP operations. Each sweep has
12 routing layers; successive sweeps traverse the layers in opposite directions.
The X- and Z-memory experiments use the same transport schedule with the
appropriate boundary preparation and logical observables.

The default experiment contains four noisy syndrome-extraction sweeps between
two ideal boundary sweeps. It uses the following circuit-level noise model,
with physical error parameter `p`:

- One `DEPOLARIZE2(p)` follows each complete CXSWAP operation. Deleted-check
  locations retain a SWAP and receive the same channel.
- Auxiliary preparation and readout use `DEPOLARIZE1(p)`. Basis rotations are
  ideal.
- Each data qubit receives one `DEPOLARIZE1(p)` idle channel per noisy sweep,
  representing the combined measurement and preparation interval.

The constructor checks a small noiseless sample and records a SHA-256 hash of
the canonical Stim circuit text. This is an implementation check, not a code-
distance computation or a statistical performance estimate.

Canonical reference hashes for both memory bases of the three highlighted codes
are supplied in [`circuit_references.json`](../src/armillary/data/circuit_references.json).
They use four noisy rounds and `p = 0.001`. A circuit hash is computed from
`str(circuit).encode("utf-8")`; file checksums are recorded separately because
the saved `.stim` file includes a final newline.

## Deleted-check routing

At deleted-check locations, the simulator applies SWAP followed by
`DEPOLARIZE2(p)`. The manuscript's zero-state explanation is an ideal-gate
identity: if the routing ancilla is the CX control and enters in an error-free
`|0>` state, CX acts trivially on any data state. Thus CXSWAP and SWAP give the
same output after the same post-gate depolarizing channel. The marginal channel
on the transported data is `DEPOLARIZE1(4p/5)` in Stim's convention.

This identity assumes the zero-state condition at each interaction. A preparation
error or an error carried by a reused ancilla can propagate through a subsequent
CX, so equal post-gate depolarizing channels alone do not make the full noisy
implementations equivalent. The reported simulations use the SWAP-only primitive
at deleted-check locations, with the preparation and gate noise specified above.

## Detector-error representation

Stim converts the circuit to a detector error model without graphlike
decomposition. Each error instruction becomes one binary fault variable with
its associated probability. `H` maps a fault pattern to detector parities;
`logical_map` records the corresponding logical-observable parities. An error
instruction containing multiple components remains a single variable. Faults
with no detector support are retained, and faults with the same syndrome but
different logical effects remain distinct.

The decoder receives only `H`, the probabilities, and a measured syndrome.
The logical map and sampled observable outcomes are used separately to score
the correction. Maximizing the probability of an individual fault pattern in
this representation is distinct from summing over all patterns in a logical
equivalence class.

## Saved model format

`DetectorModel.save(directory)` writes:

| File | Contents |
|---|---|
| `model.stim` | Physical circuit in Stim text format |
| `H.npz` | Sparse binary matrix, detectors × fault variables |
| `probabilities.npy` | One probability per fault variable |
| `logical_map.npy` | One packed logical-effect row per fault variable |
| `metadata.json` | Code and circuit parameters, dimensions and file hashes |

Logical effects are stored as `uint8` bytes with little-endian bit order; unused
bits in the last byte are zero. Saved arrays are loaded without Python pickles.
`load_model(directory)` verifies checksums and dimensions. Model files can be
reused across sampling and decoding without reconstructing the circuit.

```python
from armillary.models import build_model

model = build_model("8749", "x", rounds=4, p=0.001)
model.save("outputs/8749-x-model")
detectors, observables = model.circuit.compile_detector_sampler(seed=2026).sample(
    10, separate_observables=True
)
```

Sampling uses the physical circuit, rather than substituting independent draws
from the decoder's detector error model. Random seeds and software versions
should be recorded when comparing runs.
