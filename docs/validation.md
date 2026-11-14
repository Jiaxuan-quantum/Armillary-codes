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

The local suite passes **35 tests**, including 122 candidate subchecks. Decoder
validation used Python 3.12, NumPy 2.2.6, SciPy 1.14.1, Relay-BP 0.2.2, and COPT
8.0.7. Tests involving optional dependencies are skipped when those dependencies
are absent. The GitHub workflow tests the portable installation on Linux.

The circuit reference check uses the archived four-round circuits at p=0.001.
These checks compare canonical Stim text and detector-model dimensions, together
with noiseless consistency and sampled map identities. They are construction
checks, not new Monte Carlo performance measurements.

Distance searches and the full benchmark sampling campaign are not rerun by the
test suite. A short smoke experiment verifies that the public sampling interface
runs; it does not independently measure the rare failure rates in the benchmark.
