"""Exercise artifact identity and two-stage accounting on a one-fault circuit."""
import json
from pathlib import Path

import numpy as np
import pytest
import stim

from armillary import experiments
from armillary.models import DetectorModel


class HandoffNonzero:
    """A test front stage that leaves nonzero syndromes to the real MLE solver."""

    def __init__(self, h, probabilities):
        self.faults = h.shape[1]

    def decode(self, syndrome):
        accepted = not np.any(syndrome)
        return {"accepted": accepted,
                "correction": np.zeros(self.faults, np.uint8) if accepted else None}


@pytest.fixture
def tiny_model(tmp_path):
    # The measured error has both a detector and a logical effect. MLE has a
    # unique feasible correction for each syndrome, making scoring unambiguous.
    circuit = stim.Circuit("""
        R 0
        X_ERROR(0.25) 0
        M 0
        DETECTOR rec[-1]
        OBSERVABLE_INCLUDE(0) rec[-1]
    """)
    model = DetectorModel.from_circuit(circuit, dict(
        candidate_id="test_one_fault", basis="z", p=0.25, rounds=1, k=1))
    directory = tmp_path / "model"
    model.save(directory)
    return directory


@pytest.fixture
def front_run(tiny_model, tmp_path, monkeypatch):
    monkeypatch.setattr(experiments, "RelayDecoder", HandoffNonzero)
    output = tmp_path / "front"
    info = experiments.sample(tiny_model, 64, 72105, output, batch_size=16)
    assert 0 < info["counts"]["abstain"] < 64
    return output, info


def read_summary(directory):
    return json.loads((Path(directory) / "summary.json").read_text(encoding="utf-8"))


def test_complete_two_stage_accounting(tiny_model, front_run, tmp_path):
    directory, front = front_run
    counts = front["counts"]
    assert counts["completed"] == counts["correct"] + counts["abstain"] == 64
    assert counts["wrong"] == counts["pending"] == 0
    with np.load(directory / "handoffs.npz", allow_pickle=False) as data:
        assert str(data["run_id"].item()) == front["run_id"]
        assert len(data["sample_ids"]) == counts["abstain"]
        assert "observed_logicals" not in data.files
    result = experiments.decode_handoffs(tiny_model, directory, tmp_path / "mle")
    assert result["planned"] == result["completed"] == result["correct"] == counts["abstain"]
    assert result["wrong"] == result["unresolved"] == result["pending"] == 0
    assert result["combined"]["samples"] == 64
    assert result["combined"]["failures"] == 0
    assert result["combined"]["observed_bit_error_contribution"] == 0
    assert result["combined"]["zero_event_one_sided_95_upper"] > 0


@pytest.mark.parametrize("name", ["handoffs.npz", "scoring.npz"])
def test_modified_checkpoint_rejected(tiny_model, front_run, tmp_path, name):
    directory, _ = front_run
    with (directory / name).open("ab") as stream:
        stream.write(b"modified")
    with pytest.raises(ValueError, match="checksum"):
        experiments.decode_handoffs(tiny_model, directory, tmp_path / "mle")


def test_different_run_identity_rejected(tiny_model, front_run, tmp_path):
    directory, front = front_run
    # Even identical model/sample ordinals must not permit cross-run scoring.
    front["run_id"] = "another-run"
    (directory / "summary.json").write_text(json.dumps(front), encoding="utf-8")
    with pytest.raises(ValueError, match="identifier mismatch"):
        experiments.decode_handoffs(tiny_model, directory, tmp_path / "mle")


def test_inconsistent_counts_rejected(tiny_model, front_run, tmp_path):
    directory, front = front_run
    front["counts"]["correct"] += 1
    (directory / "summary.json").write_text(json.dumps(front), encoding="utf-8")
    with pytest.raises(ValueError, match="counts"):
        experiments.decode_handoffs(tiny_model, directory, tmp_path / "mle")


@pytest.mark.parametrize("error_type,state", [(RuntimeError, "error"),
                                             (KeyboardInterrupt, "interrupted")])
def test_front_interruption_retains_pending(tiny_model, tmp_path, monkeypatch, error_type, state):
    class StopAfterFour(HandoffNonzero):
        calls = 0

        def decode(self, syndrome):
            self.calls += 1
            if self.calls == 5:
                raise error_type("test interruption")
            return super().decode(syndrome)

    monkeypatch.setattr(experiments, "RelayDecoder", StopAfterFour)
    directory = tmp_path / "interrupted_front"
    with pytest.raises(error_type, match="test interruption"):
        experiments.sample(tiny_model, 16, 72105, directory, batch_size=16)
    front = read_summary(directory)
    assert front["state"] == state
    counts = front["counts"]
    assert counts["completed"] == counts["correct"] + counts["wrong"] + counts["abstain"] == 4
    assert counts["pending"] == 12
    result = experiments.decode_handoffs(tiny_model, directory, tmp_path / "mle")
    assert result["completed"] == counts["abstain"]
    assert "combined" not in result


def test_unresolved_mle_is_not_success(tiny_model, front_run, tmp_path, monkeypatch):
    class NoIncumbent:
        def __init__(self, *args, **kwargs):
            pass

        def decode(self, syndrome):
            return {"accepted": False, "correction": None,
                    "diagnostics": {"status": "interrupted", "optimal": False}}

    monkeypatch.setattr(experiments, "MLEDecoder", NoIncumbent)
    directory, front = front_run
    result = experiments.decode_handoffs(tiny_model, directory, tmp_path / "mle")
    assert result["unresolved"] == result["completed"] == front["counts"]["abstain"]
    assert result["correct"] == result["wrong"] == 0
    assert "combined" not in result


def test_mle_execution_error_retains_pending(tiny_model, front_run, tmp_path, monkeypatch):
    class BrokenMLE:
        def __init__(self, *args, **kwargs):
            pass

        def decode(self, syndrome):
            raise RuntimeError("test optimizer error")

    monkeypatch.setattr(experiments, "MLEDecoder", BrokenMLE)
    directory, front = front_run
    output = tmp_path / "mle"
    with pytest.raises(RuntimeError, match="test optimizer error"):
        experiments.decode_handoffs(tiny_model, directory, output)
    result = read_summary(output)
    assert result["state"] == "error"
    assert result["completed"] == result["correct"] == result["wrong"] == result["unresolved"] == 0
    assert result["pending"] == front["counts"]["abstain"]
    assert "combined" not in result
