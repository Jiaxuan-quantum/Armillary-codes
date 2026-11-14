"""Statistical conventions and reconciliation of the released table."""
from importlib.resources import files
import json
import math

import pytest

from armillary.statistics import equivalent_rate, failure_statistics


def test_zero_events_have_positive_upper_bound():
    result = failure_statistics(0, 1000, 100, 4)
    assert result["equivalent_rate"] == 0
    assert result["zero_event_one_sided_95_upper"] == pytest.approx(1 - .05 ** (1 / 1000))
    assert result["equivalent_zero_event_one_sided_95_upper"] > 0


def test_small_equivalent_rate_keeps_precision():
    rate = equivalent_rate(1 / 35_484_090, 4348, 4)
    assert rate == pytest.approx(1.6203797505343838e-12, rel=1e-12, abs=0)
    assert -math.expm1(4348 * 4 * math.log1p(-rate)) == pytest.approx(1 / 35_484_090)


@pytest.mark.parametrize("failures,samples", [(-1, 10), (11, 10), (0, 0)])
def test_invalid_counts_rejected(failures, samples):
    with pytest.raises(ValueError):
        failure_statistics(failures, samples, 1, 4)


def test_released_benchmark_reconciles():
    source = files("armillary").joinpath("data", "benchmark.json")
    data = json.loads(source.read_text(encoding="utf-8"))
    assert len(data["pooled"]) == 3 and len(data["sides"]) == 6
    sum_fields = ("N", "C1", "W1", "B", "C2", "W2", "bits1", "bits2")
    expected_failures = {"8749": 24, "9053": 10, "p181": 1}
    for pooled in data["pooled"]:
        sides = [r for r in data["sides"] if r["candidate"] == pooled["candidate"]]
        assert {r["basis"] for r in sides} == {"x", "z"}
        assert sides[0]["N"] == sides[1]["N"]
        for field in sum_fields:
            assert pooled[field] == sum(r[field] for r in sides)
        for row in [pooled, *sides]:
            assert row["N"] == row["C1"] + row["W1"] + row["B"]
            assert row["B"] == row["C2"] + row["W2"]
        assert pooled["W1"] + pooled["W2"] == expected_failures[pooled["candidate"]]
