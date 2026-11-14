"""Detector-error models and their syndrome/logical maps.

Each Stim error instruction is retained as one binary fault variable. In
particular, columns with the same syndrome but different logical effects are
distinct, and faults with no detector support are not discarded.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csc_matrix, load_npz, save_npz
import stim

from .candidates import load_candidate
from .circuits import build_circuit

SCHEMA = "armillary.detector_model.v1"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _dem_arrays(dem, detector_count, observable_count):
    """Preserve the flattened DEM's column order and complete fault effects."""
    dem = dem.flattened()
    logical = np.zeros((dem.num_errors, (observable_count + 7) // 8), np.uint8)
    indices, indptr, priors = [], [0], []
    column = 0
    for instruction in dem:
        if instruction.type != "error":
            continue
        detectors, observables = set(), set()
        for target in instruction.targets_copy():
            if target.is_relative_detector_id():
                detectors.symmetric_difference_update((target.val,))
            elif target.is_logical_observable_id():
                observables.symmetric_difference_update((target.val,))
            # A separator groups components of the same fault. It does not
            # create independent variables; the full instruction is one column.
            elif not target.is_separator():
                raise ValueError("Unsupported detector-model target")
        if any(v >= detector_count for v in detectors):
            raise ValueError("Detector target outside the circuit")
        if any(v >= observable_count for v in observables):
            raise ValueError("Observable target outside the circuit")
        indices.extend(sorted(detectors))
        indptr.append(len(indices))
        priors.append(instruction.args_copy()[0])
        for observable in observables:
            logical[column, observable // 8] ^= np.uint8(1 << (observable % 8))
        column += 1
    h = csc_matrix(
        (np.ones(len(indices), np.uint8), np.asarray(indices, np.int32),
         np.asarray(indptr, np.int64)), shape=(detector_count, column))
    return h, np.asarray(priors, np.float64), logical


@dataclass
class DetectorModel:
    """Circuit, binary detector matrix, priors, and packed logical effects.

    ``H`` has shape (detectors, fault variables). ``logical_map`` stores one
    little-endian packed observable vector per fault variable. Logical truth
    from a sampled experiment is never an input to a decoder.
    """

    circuit: stim.Circuit
    H: csc_matrix
    probabilities: np.ndarray
    logical_map: np.ndarray
    metadata: dict

    def __post_init__(self):
        self.H = csc_matrix(self.H, copy=False)
        self.H.sum_duplicates()
        self.H.eliminate_zeros()
        self.H.sort_indices()
        if np.any(self.H.data != 1):
            raise ValueError("The detector matrix must be binary")
        self.H = self.H.astype(np.uint8, copy=False)
        self.probabilities = np.asarray(self.probabilities, dtype=np.float64)
        faults = self.H.shape[1]
        observables = self.circuit.num_observables
        if self.H.shape[0] != self.circuit.num_detectors:
            raise ValueError("Detector matrix and circuit dimensions differ")
        if (self.probabilities.shape != (faults,)
                or not np.all(np.isfinite(self.probabilities))
                or np.any(self.probabilities <= 0)
                or np.any(self.probabilities >= 0.5)):
            raise ValueError("Each fault prior must be finite and in (0, 0.5)")
        if (self.logical_map.dtype != np.uint8
                or self.logical_map.shape != (faults, (observables + 7) // 8)):
            raise ValueError("Invalid packed logical map")
        if observables % 8 and np.any(self.logical_map[:, -1] >> (observables % 8)):
            raise ValueError("Logical-map padding bits must be zero")

    def predict_observables(self, correction) -> np.ndarray:
        """Map a binary fault pattern to its logical-observable flips."""
        correction = np.asarray(correction)
        if (correction.shape != (self.H.shape[1],)
                or not np.all((correction == 0) | (correction == 1))):
            raise ValueError("Correction must be a binary vector with one entry per fault")
        selected = np.flatnonzero(correction)
        packed = (np.bitwise_xor.reduce(self.logical_map[selected], axis=0)
                  if len(selected) else np.zeros(self.logical_map.shape[1], np.uint8))
        return np.unpackbits(packed, bitorder="little")[:self.circuit.num_observables]

    def save(self, directory) -> Path:
        """Save portable arrays and file checksums for subsequent decoding."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        metadata_path = directory / "metadata.json"
        if metadata_path.exists():
            previous = json.loads(metadata_path.read_text(encoding="utf-8"))
            if previous.get("model_id") != self.metadata.get("model_id"):
                raise FileExistsError("The directory already contains a different model")
            # Avoid rewriting arrays that may be memory-mapped by this object.
            expected = previous.get("files", {})
            if expected and all((directory / n).is_file() and _sha(directory / n) == h
                                for n, h in expected.items()):
                return directory
            raise FileExistsError("The existing model is incomplete or modified; use a new directory")
        # Python handles Unicode paths consistently on Windows as well as POSIX.
        (directory / "model.stim").write_text(str(self.circuit) + "\n", encoding="utf-8", newline="\n")
        save_npz(directory / "H.npz", self.H)
        np.save(directory / "probabilities.npy", self.probabilities, allow_pickle=False)
        np.save(directory / "logical_map.npy", self.logical_map, allow_pickle=False)
        metadata = dict(self.metadata)
        metadata["files"] = {name: _sha(directory / name) for name in
                             ("model.stim", "H.npz", "probabilities.npy", "logical_map.npy")}
        metadata_path.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n",
                                 encoding="utf-8", newline="\n")
        return directory

    @classmethod
    def from_circuit(cls, circuit, metadata=None):
        """Translate a circuit using Stim without graphlike decomposition."""
        dem = circuit.detector_error_model(
            decompose_errors=False, allow_gauge_detectors=False).flattened()
        h, priors, logical = _dem_arrays(dem, circuit.num_detectors, circuit.num_observables)
        info = dict(metadata or {})
        info.update(schema=SCHEMA, stim_version=stim.__version__,
                    detector_count=circuit.num_detectors,
                    observable_count=circuit.num_observables, fault_variables=h.shape[1],
                    circuit_sha256=hashlib.sha256(str(circuit).encode()).hexdigest(),
                    dem_sha256=hashlib.sha256(str(dem).encode()).hexdigest())
        info["model_id"] = hashlib.sha256(
            (SCHEMA + "\n" + info["circuit_sha256"] + "\n" + info["dem_sha256"]).encode()
        ).hexdigest()
        model = cls(circuit, h, priors, logical, info)
        # This fixed tiny check validates both maps against Stim's own sampler.
        detector, observable, fault = dem.compile_sampler(seed=26110401).sample(
            8, return_errors=True)
        for d, o, e in zip(detector, observable, fault):
            if not np.array_equal(np.asarray(h @ e.astype(np.uint8)).ravel() & 1, d):
                raise ValueError("Detector-model syndrome mapping failed")
            if not np.array_equal(model.predict_observables(e), o):
                raise ValueError("Detector-model logical mapping failed")
        return model


def build_model(identifier, basis, rounds=4, p=0.001) -> DetectorModel:
    """Reconstruct a catalog candidate and its memory experiment."""
    record = load_candidate(identifier)
    circuit, audit = build_circuit(record, basis, rounds=rounds, p=p)
    if not audit["noiseless_pass"]:
        raise ValueError("The noiseless memory experiment is inconsistent")
    metadata = dict(candidate_id=record["candidate_id"], n=record["n"], k=record["k"],
                    distance=record.get("exact_distance"), basis=basis, rounds=rounds,
                    p=float(p), matrix_sha256=record["matrix_sha256"], circuit_audit=audit,
                    noise_model="depolarizing-preparation-readout-single-data-idle")
    return DetectorModel.from_circuit(circuit, metadata)


def load_model(directory, verify=True) -> DetectorModel:
    """Read a saved model, checking its files before exposing the decoder inputs."""
    directory = Path(directory)
    info = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    if info.get("schema") != SCHEMA:
        raise ValueError("Unsupported detector-model schema")
    expected_names = {"model.stim", "H.npz", "probabilities.npy", "logical_map.npy"}
    if set(info.get("files", {})) != expected_names:
        raise ValueError("Incomplete detector-model file manifest")
    if verify:
        for name, expected in info["files"].items():
            if _sha(directory / name) != expected:
                raise ValueError("Detector-model checksum mismatch: " + name)
    model = DetectorModel(
        stim.Circuit((directory / "model.stim").read_text(encoding="utf-8")),
        load_npz(directory / "H.npz"),
        np.load(directory / "probabilities.npy", allow_pickle=False),
        np.load(directory / "logical_map.npy", mmap_mode="r", allow_pickle=False), info)
    if (model.H.shape != (info["detector_count"], info["fault_variables"])
            or model.circuit.num_observables != info["observable_count"]):
        raise ValueError("Model dimensions disagree with metadata")
    return model
