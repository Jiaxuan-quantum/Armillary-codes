"""Small, portable memory experiments and saved handoffs to MLE.

Decoder inputs contain only syndromes and the fixed detector model. Simulated
logical outcomes are stored separately and used only after a correction exists.
"""
from __future__ import annotations

import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import platform
import uuid

import numpy as np

from .decoders import MLEDecoder, RelayDecoder
from .models import load_model
from .statistics import failure_statistics


def _new_directory(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError("Choose an empty output directory: " + str(path))
    path.mkdir(parents=True, exist_ok=True)
    return path


def _json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _npz(path, **arrays):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    temporary.replace(path)


def _sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _runtime(*packages):
    versions = {"python": platform.python_version()}
    for package in packages:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    return versions


def sample(model_directory, shots, seed, output, batch_size=32):
    """Sample a fixed number of circuit shots and apply the Relay front stage.

    Completed outcomes are checkpointed after each batch. If interrupted, the
    summary explicitly separates completed shots from the requested remainder.
    Repeated seeds with identical circuit/runtime/batching repeat samples.
    """
    if shots < 1 or batch_size < 1 or seed < 0:
        raise ValueError("shots/batch_size must be positive and seed nonnegative")
    model = load_model(model_directory)
    decoder = RelayDecoder(model.H, model.probabilities)
    out = _new_directory(output)
    sampler = model.circuit.compile_detector_sampler(seed=seed)
    counts = dict(requested=shots, completed=0, correct=0, wrong=0, abstain=0,
                  logical_bit_errors=0, pending=shots)
    handoffs, truths, sample_ids = [], [], []
    run_id = str(uuid.uuid4())
    info = dict(schema="armillary.front_run.v1", run_id=run_id,
                model_id=model.metadata["model_id"], candidate_id=model.metadata["candidate_id"],
                basis=model.metadata["basis"], p=model.metadata["p"],
                rounds=model.metadata["rounds"], k=model.metadata["k"],
                runtime=_runtime("armillary-codes", "numpy", "stim", "relay-bp"),
                seed=seed, batch_size=batch_size, state="running", counts=counts)

    def checkpoint():
        counts["pending"] = shots - counts["completed"]
        _npz(out / "handoffs.npz", model_id=model.metadata["model_id"], run_id=run_id,
            sample_ids=np.asarray(sample_ids, dtype=np.int64),
            syndromes=np.asarray(handoffs, dtype=np.uint8).reshape(-1, model.H.shape[0]))
        _npz(out / "scoring.npz", model_id=model.metadata["model_id"], run_id=run_id,
            sample_ids=np.asarray(sample_ids, dtype=np.int64),
            observed_logicals=np.asarray(truths, dtype=np.uint8).reshape(-1, model.circuit.num_observables))
        info["files"] = {name: _sha(out / name) for name in ("handoffs.npz", "scoring.npz")}
        _json(out / "summary.json", info)

    try:
        for start in range(0, shots, batch_size):
            size = min(batch_size, shots - start)
            syndromes, observed = sampler.sample(size, separate_observables=True)
            for offset, (syndrome, truth) in enumerate(zip(syndromes, observed)):
                decoded = decoder.decode(syndrome)
                if decoded["accepted"]:
                    errors = int(np.count_nonzero(model.predict_observables(decoded["correction"]) != truth))
                    counts["wrong" if errors else "correct"] += 1
                    counts["logical_bit_errors"] += errors
                else:
                    counts["abstain"] += 1
                    handoffs.append(syndrome.copy())
                    truths.append(truth.copy())
                    sample_ids.append(start + offset)
                counts["completed"] += 1
            checkpoint()
            print(f"{counts['completed']}/{shots}: correct={counts['correct']} "
                  f"wrong={counts['wrong']} handoffs={counts['abstain']}", flush=True)
    except BaseException as error:
        info.update(state="interrupted" if isinstance(error, KeyboardInterrupt) else "error",
                    error_type=type(error).__name__)
        checkpoint()
        raise
    info["state"] = "complete"
    checkpoint()
    return info


def decode_handoffs(model_directory, front_directory, output, backend="scipy"):
    """Apply MLE to a saved handoff set, scoring corrections separately."""
    model = load_model(model_directory)
    front_directory = Path(front_directory)
    front = json.loads((front_directory / "summary.json").read_text(encoding="utf-8"))
    model_id = model.metadata["model_id"]
    if front["model_id"] != model_id:
        raise ValueError("Handoffs were produced with a different model")
    for name in ("handoffs.npz", "scoring.npz"):
        if _sha(front_directory / name) != front.get("files", {}).get(name):
            raise ValueError("Handoff checkpoint checksum mismatch: " + name)
    counts = front["counts"]
    if (counts["completed"] != counts["correct"] + counts["wrong"] + counts["abstain"]
            or counts["pending"] != counts["requested"] - counts["completed"]
            or any(not isinstance(value, int) or value < 0 for value in counts.values())
            or any(front[key] != model.metadata[key] for key in ("candidate_id", "basis", "p", "k", "rounds"))):
        raise ValueError("Inconsistent first-stage counts or model metadata")
    with np.load(front_directory / "handoffs.npz", allow_pickle=False) as data:
        if str(data["model_id"].item()) != model_id or str(data["run_id"].item()) != front["run_id"]:
            raise ValueError("Handoff model identifier mismatch")
        syndromes, ids = data["syndromes"].copy(), data["sample_ids"].copy()
    with np.load(front_directory / "scoring.npz", allow_pickle=False) as data:
        if (str(data["model_id"].item()) != model_id or str(data["run_id"].item()) != front["run_id"]
                or not np.array_equal(ids, data["sample_ids"])):
            raise ValueError("Scoring data do not match the handoff set")
        truths = data["observed_logicals"].copy()
    if (ids.ndim != 1 or not np.issubdtype(ids.dtype, np.integer)
            or np.any(ids < 0) or np.any(ids >= counts["completed"])
            or syndromes.shape != (len(ids), model.H.shape[0])
            or truths.shape != (len(ids), model.circuit.num_observables)
            or len(np.unique(ids)) != len(ids)
            or len(ids) != front["counts"]["abstain"]
            or not np.all((syndromes == 0) | (syndromes == 1))
            or not np.all((truths == 0) | (truths == 1))):
        raise ValueError("Invalid handoff/scoring arrays")
    decoder = MLEDecoder(model.H, model.probabilities, backend=backend)
    out = _new_directory(output)
    info = dict(schema="armillary.mle_run.v1", front_run_id=front["run_id"],
                model_id=model_id, backend=backend, state="running", planned=len(ids),
                runtime=_runtime("armillary-codes", "numpy", "scipy", "coptpy"),
                completed=0, correct=0, wrong=0, unresolved=0, logical_bit_errors=0,
                pending=len(ids))
    rows = []

    def checkpoint():
        info["pending"] = len(ids) - info["completed"]
        _json(out / "outcomes.json", rows)
        info["outcomes_sha256"] = _sha(out / "outcomes.json")
        _json(out / "summary.json", info)

    try:
        for sample_id, syndrome, truth in zip(ids, syndromes, truths):
            decoded = decoder.decode(syndrome)
            row = dict(sample_id=int(sample_id), accepted=bool(decoded["accepted"]),
                       diagnostics=decoded["diagnostics"])
            if decoded["accepted"]:
                correction = decoded["correction"]
                errors = int(np.count_nonzero(model.predict_observables(correction) != truth))
                info["wrong" if errors else "correct"] += 1
                info["logical_bit_errors"] += errors
                row.update(correct=errors == 0, logical_bit_errors=errors,
                           fault_indices=np.flatnonzero(correction).tolist())
            else:
                info["unresolved"] += 1
                row["correct"] = None
            rows.append(row)
            info["completed"] += 1
            checkpoint()
            print(f"MLE {info['completed']}/{len(ids)}: correct={info['correct']} "
                  f"wrong={info['wrong']} unresolved={info['unresolved']}", flush=True)
    except BaseException as error:
        info.update(state="interrupted" if isinstance(error, KeyboardInterrupt) else "error",
                    error_type=type(error).__name__)
        checkpoint()
        raise
    info["state"] = "complete"
    if not info["unresolved"] and not front["counts"]["pending"]:
        failures = front["counts"]["wrong"] + info["wrong"]
        info["combined"] = failure_statistics(failures, front["counts"]["completed"],
                                              front["k"], front["rounds"])
        info["combined"]["observed_bit_error_contribution"] = (
            (front["counts"]["logical_bit_errors"] + info["logical_bit_errors"])
            / (front["counts"]["completed"] * front["k"] * front["rounds"]))
    checkpoint()
    return info
