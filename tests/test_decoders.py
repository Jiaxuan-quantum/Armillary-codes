"""Check decoder semantics and the parity formulation on tractable examples."""
import itertools
import sys
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from armillary.decoders import MLEDecoder, RelayDecoder


def test_mle_matches_exhaustive_patterns():
    # The fourth column is detectorless and must be retained in the model.
    h = np.array([[1, 1, 0, 0], [0, 1, 1, 0]], dtype=np.uint8)
    probabilities = np.array([0.11, 0.18, 0.07, 0.05])
    decoder = MLEDecoder(h, probabilities)
    weights = np.log((1 - probabilities) / probabilities)
    for syndrome in itertools.product((0, 1), repeat=h.shape[0]):
        syndrome = np.array(syndrome, dtype=np.uint8)
        feasible = [np.array(e, dtype=np.uint8)
                    for e in itertools.product((0, 1), repeat=h.shape[1])
                    if np.array_equal((h @ np.array(e)) % 2, syndrome)]
        target = min(float(weights @ e) for e in feasible)
        result = decoder.decode(syndrome)
        assert result["accepted"]
        assert result["diagnostics"]["optimal"]
        assert result["diagnostics"]["objective"] == pytest.approx(target)
        assert np.array_equal(h @ result["correction"] % 2, syndrome)
        assert result["correction"][3] == 0


def test_mle_infeasible_is_not_success():
    result = MLEDecoder([[1], [1]], [0.1]).decode([1, 0])
    assert result["correction"] is None
    assert not result["accepted"]
    assert result["diagnostics"]["status"] == "infeasible"


@pytest.mark.parametrize("decoder_type", [RelayDecoder, MLEDecoder])
def test_bad_inputs_rejected(decoder_type):
    with pytest.raises(ValueError):
        decoder_type([[2]], [0.1])
    with pytest.raises(ValueError):
        decoder_type([[1]], [0.5])
    with pytest.raises(ValueError):
        decoder_type([[1]], [0.1]).decode([2])


def _relay_mock(outputs, calls):
    class FakeRelay:
        def __init__(self, h, error_priors, **parameters):
            calls.append(parameters)

        def decode_detailed(self, syndrome):
            return SimpleNamespace(decoding=np.array(outputs.pop(0), dtype=np.uint8),
                                   iterations=5)
    return SimpleNamespace(RelayDecoderF64=FakeRelay)


def test_relay_short_success_never_runs_long():
    calls = []
    with patch.dict(sys.modules, {"relay_bp": _relay_mock([[1, 0]], calls)}):
        result = RelayDecoder([[1, 1]], [0.1, 0.1]).decode([1])
    assert result["accepted"]
    assert result["diagnostics"]["selected_stage"] == "short"
    assert calls == [dict(pre_iter=200, num_sets=20, set_max_iter=100,
                          gamma0=0.125, stop_nconv=2, seed=0)]


def test_relay_invalid_short_runs_fresh_long():
    calls = []
    with patch.dict(sys.modules, {"relay_bp": _relay_mock([[0, 0], [0, 1]], calls)}):
        result = RelayDecoder([[1, 1]], [0.1, 0.1]).decode([1])
    assert result["accepted"]
    assert result["diagnostics"]["selected_stage"] == "long"
    assert calls[1] == dict(pre_iter=500, num_sets=200, set_max_iter=200,
                            gamma0=0.35, stop_nconv=2, seed=0)


def test_relay_abstains_if_both_outputs_invalid():
    with patch.dict(sys.modules, {"relay_bp": _relay_mock([[0], [0]], [])}):
        result = RelayDecoder([[1]], [0.1]).decode([1])
    assert not result["accepted"]
    assert result["correction"] is None


def test_relay_zero_syndrome_is_explicit_fast_path():
    with patch.dict(sys.modules, {"relay_bp": None}):
        result = RelayDecoder([[1]], [0.1]).decode([0])
    assert result["accepted"]
    assert result["diagnostics"]["selected_stage"] == "zero_syndrome"


def test_relay_library_exception_is_not_abstention():
    class BrokenRelay:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("native failure")
    with patch.dict(sys.modules, {"relay_bp": SimpleNamespace(RelayDecoderF64=BrokenRelay)}):
        with pytest.raises(RuntimeError, match="native failure"):
            RelayDecoder([[1]], [0.1]).decode([1])


def test_real_relay_repetition_example():
    pytest.importorskip("relay_bp")
    h = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.uint8)
    result = RelayDecoder(h, [0.1, 0.1, 0.1]).decode([1, 1])
    assert result["accepted"]
    assert result["correction"].tolist() == [0, 1, 0]


def test_copt_matches_exhaustive_patterns():
    pytest.importorskip("coptpy")
    h = np.array([[1, 1, 0, 0], [0, 1, 1, 0]], dtype=np.uint8)
    probabilities = np.array([0.11, 0.18, 0.07, 0.05])
    decoder = MLEDecoder(h, probabilities, backend="copt")
    weights = np.log((1 - probabilities) / probabilities)
    for s in itertools.product((0, 1), repeat=2):
        syndrome = np.array(s, dtype=np.uint8)
        feasible = [np.array(e, dtype=np.uint8)
                    for e in itertools.product((0, 1), repeat=4)
                    if np.array_equal(h @ np.array(e) % 2, syndrome)]
        result = decoder.decode(syndrome)
        assert result["accepted"]
        assert result["diagnostics"]["optimal"]
        assert result["diagnostics"]["objective"] == pytest.approx(
            min(float(weights @ e) for e in feasible))
    result = MLEDecoder([[1], [1]], [0.1], backend="copt").decode([1, 0])
    assert not result["accepted"]
    assert result["diagnostics"]["status"] == "infeasible"
