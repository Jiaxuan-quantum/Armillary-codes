"""Block-failure statistics with an explicit reporting convention."""
from __future__ import annotations

import math

from scipy.stats import beta


def equivalent_rate(block_rate: float, k: int, rounds: int) -> float:
    """Block-normalized equivalent per-logical-qubit, per-round probability."""
    if not 0 <= block_rate <= 1 or k <= 0 or rounds <= 0:
        raise ValueError("Expected a probability and positive k and rounds")
    return 1.0 if block_rate == 1 else -math.expm1(math.log1p(-block_rate) / (k * rounds))


def failure_statistics(failures: int, samples: int, k: int, rounds: int) -> dict:
    """Return estimates and nominal fixed-N binomial reference intervals.

    These intervals do not adjust for adaptive stopping or pooled heterogeneous
    memory experiments. Zero observed failures do not establish zero risk.
    """
    if samples <= 0 or not 0 <= failures <= samples:
        raise ValueError("Expected 0 <= failures <= samples with samples > 0")
    lo = float(beta.ppf(.025, failures, samples - failures + 1)) if failures else 0.0
    hi = float(beta.ppf(.975, failures + 1, samples - failures)) if failures < samples else 1.0
    rate = failures / samples
    result = dict(samples=samples, failures=failures, block_failure_rate=rate,
                  equivalent_rate=equivalent_rate(rate, k, rounds),
                  nominal_binomial_95=[lo, hi],
                  equivalent_nominal_95=[equivalent_rate(lo, k, rounds),
                                         equivalent_rate(hi, k, rounds)])
    if failures == 0:
        upper = -math.expm1(math.log(.05) / samples)
        result["zero_event_one_sided_95_upper"] = upper
        result["equivalent_zero_event_one_sided_95_upper"] = equivalent_rate(upper, k, rounds)
    return result
