# Benchmark counts and statistical conventions

`src/armillary/data/benchmark.json` contains aggregate counts from the three
representative codes at physical error probability 0.001, with four noisy rounds
and ideal boundaries. Each code has equal numbers of completed X- and Z-memory
samples. The pooled population counts each memory experiment once; it is not a
joint measurement of X and Z logical failure on the same shot.

| Code | Basis | N | C1 | W1 | B | C2 | W2 |
|---|---|---:|---:|---:|---:|---:|---:|
| 8749 | X | 3,361,646 | 3,361,385 | 13 | 248 | 248 | 0 |
| 8749 | Z | 3,361,646 | 3,361,385 | 9 | 252 | 252 | 0 |
| 9053 | X | 10,386,959 | 10,386,478 | 3 | 478 | 478 | 0 |
| 9053 | Z | 10,386,959 | 10,386,431 | 6 | 522 | 522 | 0 |
| p181 | X | 17,742,045 | 17,740,645 | 0 | 1,400 | 1,400 | 0 |
| p181 | Z | 17,742,045 | 17,740,002 | 1 | 2,042 | 2,042 | 0 |

- `N`: completed first-stage experiments.
- `C1`, `W1`: correct and incorrect accepted Relay-BP outputs.
- `B`: Relay abstentions handed to MLE.
- `C2`, `W2`: correct and incorrect MLE outputs.

Here `N = C1 + W1 + B`, `B = C2 + W2`, and the final block-failure count is
`F = W1 + W2`. Syndrome consistency does not imply correct logical action. An
incorrect accepted Relay output remains a failure even though it never enters
MLE. An execution error or an experiment without a completed output is neither a
successful decode nor a Relay abstention.

The p181 population combines several sampling batches. Equal X/Z exposure excludes
192 completed but unpaired samples, all correctly decoded by Relay. Another 192
samples without a completed first-stage result are excluded. Neither exclusion
changes the failure or handoff counts. All 4,942 included handoffs have a scored
final correction. The MLE columns group outcomes by correctness; numerical global
optimality was not certified for every historical correction.

## Reported probabilities

The empirical block failure probability is `P_fail = F/N`. For a code with `k`
logical qubits and `r` noisy rounds, the manuscript reports

```text
P_L = 1 - (1 - P_fail)^(1/(k r)).
```

This equivalent probability normalizes a block event. It differs from the directly
observed logical-bit contribution, `sum(logical bit flips)/(k r N)`, because one
block failure can affect several logical observables. The latter values are
3.92e-7, 1.87e-8, and 1.46e-11 for 8749, 9053, and p181, respectively.

`armillary benchmark` recomputes point estimates and two-sided 95% Clopper–Pearson
reference intervals. It transforms both endpoints using the same equation above.
The statistics module also reports a one-sided 95% upper bound for zero events.

These are nominal fixed-size binomial intervals. Batches stopped at prescribed
handoff counts or elapsed-time cutoffs, and additional p181 data were collected
after earlier results were inspected. The pooled X/Z populations can also have
different failure probabilities. Consequently, the intervals describe counting
uncertainty; they are not a stopping-rule-adjusted coverage guarantee. The single
observed p181 failure does not certify a rate below 1e-12.

## Reproduction scope

The supplied counts reproduce the benchmark arithmetic and plot. Candidate
definitions, physical circuits, and the current Relay parameters are separately
reconstructible. The MLE implementation exposes the same error-pattern objective
as a portable reference. These commands do not replay the historical compute
campaign or guarantee the same individual solver corrections.
