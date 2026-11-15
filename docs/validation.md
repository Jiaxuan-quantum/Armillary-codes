# Validation scope

The reference release was checked with the following independent layers:

1. **Candidate reconstruction:** all 122 definitions reproduce their recorded
   ordered X/Z matrix hashes, dimensions, CSS commutation, and forward/reverse
   transport checks. All matrix identifiers are distinct.
2. **Detector-model mapping:** tests retain detectorless faults, distinguish
   identical syndromes with different logical actions, check both maps against
   Stim, and reject modified saved models.
3. **Decoders:** small parity problems are compared with exhaustive enumeration
   using both HiGHS and COPT. Actual Relay-BP 0.2.2 is exercised on a small
   syndrome. Separate dispatch tests cover the short/long cascade and abstention.
4. **Experiment accounting:** a small Stim experiment exercises saved handoffs,
   numerical MLE, and independent logical scoring. Checks cover cross-run input
   mismatches, incomplete work, execution errors, and unresolved solver outputs.
5. **Statistics and packaging:** per-basis and pooled benchmark counts reconcile;
   very small probabilities and zero-event bounds are checked. A built Python
   wheel contains the 122 definitions and can reconstruct a highlighted code
   without access to the original research directory.

The original reference-release suite passed **35 tests**, including 122 candidate subchecks. Decoder
validation used Python 3.12, NumPy 2.2.6, SciPy 1.14.1, Relay-BP 0.2.2, and COPT
8.0.7. Tests involving optional dependencies are skipped when those dependencies
are absent. The original Linux GitHub workflow also passed: 34 tests and 122 candidate
subchecks, with the optional COPT test skipped because COPT is not installed
there. COPT was tested locally as described above.

All six representative models (8749, 9053, and p181, each in both memory bases)
match the archived four-round circuits at p=0.001. Canonical circuit SHA-256,
detector-model dimensions, noiseless consistency, and sampled detector/logical
map identities pass. This comparison used Stim 1.16.0; the archived references
used Stim 1.15.0. The references are included in
[`circuit_references.json`](../src/armillary/data/circuit_references.json).
These are construction checks, not new Monte Carlo performance measurements.

Distance searches and the full benchmark sampling campaign are not rerun by the
test suite. A short smoke experiment verifies that the public sampling interface
runs; it does not independently measure the rare failure rates in the benchmark.

## Hardware-layout reproduction

The standalone `hardware_layouts` module verifies the three representative layouts separately from the circuit/decoder suite. It checks the complete fixed physical connection set, code-matrix identity, explicit continuous curves, common-endpoint exemptions, nonendpoint device regions, endpoint degree limits, and tier-averaged hardware metrics. It enumerates all stackings for each frozen face assignment. A separate GitHub job runs this reproduction on Python 3.11 and 3.14.

The archived lower-bound distance evidence is not rerun. Successful hardware verification means that the selected inputs reproduce the stated geometric bounds and resource metrics. It is not a fabrication or electromagnetic validation, a global layout-optimality proof, or a new QEC performance measurement.

For the hardware-layout addition, the local Python 3.14.4 suite passed 37 tests and 122 candidate subchecks, with two optional-dependency tests skipped. All three layouts also passed the standalone reproduction and a separate relocated-copy run. Input modification, disabled assertions, and attempts to write into the frozen input directory are rejected. These additional checks do not replace the original optional-decoder validation described above.
