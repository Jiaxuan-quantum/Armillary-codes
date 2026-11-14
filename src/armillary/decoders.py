"""Syndrome-only Relay-BP and error-pattern MLE decoders.

Logical observables and sampled logical outcomes are deliberately absent from
these interfaces. A syndrome-consistent correction need not be logically right.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import time

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, eye, hstack


@dataclass(frozen=True)
class RelayStage:
    name: str
    pre_iter: int
    num_sets: int
    set_max_iter: int
    gamma0: float
    stop_nconv: int = 2


RELAY_STAGES = (
    RelayStage("short", 200, 20, 100, 0.125),
    RelayStage("long", 500, 200, 200, 0.35),
)


def _channel(H, probabilities):
    h = csr_matrix(H, copy=True)
    h.sum_duplicates()
    h.eliminate_zeros()
    h.sort_indices()
    if h.ndim != 2 or min(h.shape) < 1:
        raise ValueError("H must have positive row and column counts")
    if not np.all(h.data == 1):
        raise ValueError("H must be binary, without duplicate nonzero entries")
    h = h.astype(np.uint8)
    probabilities = np.asarray(probabilities, dtype=np.float64).copy()
    if (probabilities.shape != (h.shape[1],)
            or not np.all(np.isfinite(probabilities))
            or not np.all((probabilities > 0) & (probabilities < 0.5))):
        raise ValueError("Each column requires one finite probability in (0, 0.5)")
    weights = np.log1p(-probabilities) - np.log(probabilities)
    return h, probabilities, weights


def _binary(value, size, label):
    value = np.asarray(value)
    if value.shape != (size,) or not np.all((value == 0) | (value == 1)):
        raise ValueError(f"{label} must be a binary vector of length {size}")
    return value.astype(np.uint8, copy=True)


def _valid(h, correction, syndrome):
    return np.array_equal(np.asarray(h @ correction).ravel() & 1, syndrome)


class RelayDecoder:
    """Short Relay followed by long Relay only if the first output is invalid.

    Uses relay-bp 0.2.2's float64 implementation. Each call constructs fresh
    decoders with seed 0, so one sample never inherits another sample's state.
    """

    def __init__(self, H, probabilities):
        self.h, self.probabilities, self.weights = _channel(H, probabilities)

    def decode(self, syndrome):
        syndrome = _binary(syndrome, self.h.shape[0], "Syndrome")
        started = time.perf_counter()
        diagnostics = {"decoder": "short_long_relay", "stages": [],
                       "selected_stage": None}
        correction = None
        if not syndrome.any():
            correction = np.zeros(self.h.shape[1], dtype=np.uint8)
            diagnostics["selected_stage"] = "zero_syndrome"
        else:
            try:
                import relay_bp
            except ImportError as exc:
                raise ImportError(
                    "Relay decoding needs relay-bp==0.2.2; install the relay extra."
                ) from exc
            for stage in RELAY_STAGES:
                parameters = asdict(stage)
                parameters.pop("name")
                decoder = relay_bp.RelayDecoderF64(
                    self.h, error_priors=self.probabilities, seed=0, **parameters)
                result = decoder.decode_detailed(syndrome)
                candidate = _binary(result.decoding, self.h.shape[1], "Correction")
                valid = _valid(self.h, candidate, syndrome)
                diagnostics["stages"].append({
                    "stage": stage.name, "syndrome_valid": bool(valid),
                    "iterations": int(result.iterations),
                })
                if valid:
                    correction = candidate
                    diagnostics["selected_stage"] = stage.name
                    break
        accepted = correction is not None
        diagnostics.update(
            reason="syndrome_valid" if accepted else "abstain",
            objective=float(self.weights @ correction) if accepted else None,
            seconds=time.perf_counter() - started,
        )
        return {"correction": correction, "accepted": accepted,
                "diagnostics": diagnostics}


class MLEDecoder:
    r"""Minimize sum_j log((1-p_j)/p_j) e_j subject to H e = s (mod 2).

    The formulation maximizes the probability of a single independent-DEM
    error pattern. It does not sum over errors in a logical equivalence class.
    SciPy/HiGHS is the portable backend; COPT is an optional licensed backend.
    """

    def __init__(self, H, probabilities, backend="scipy"):
        self.h, self.probabilities, self.weights = _channel(H, probabilities)
        if backend not in ("scipy", "copt"):
            raise ValueError("Backend must be 'scipy' or 'copt'")
        self.backend = backend
        self.row_weights = np.asarray(self.h.sum(axis=1), dtype=np.int64).ravel()

    def _scipy(self, syndrome):
        from scipy.optimize import Bounds, LinearConstraint, milp

        m, n = self.h.shape
        matrix = csc_matrix(hstack((self.h, -2 * eye(m, format="csc"))))
        # HiGHS' sparse interface expects 32-bit index arrays.
        matrix.indices = matrix.indices.astype(np.int32)
        matrix.indptr = matrix.indptr.astype(np.int32)
        upper = np.concatenate((np.ones(n), np.floor(self.row_weights / 2)))
        result = milp(
            c=np.concatenate((self.weights, np.zeros(m))),
            integrality=np.ones(n + m, dtype=np.uint8),
            bounds=Bounds(np.zeros(n + m), upper),
            constraints=LinearConstraint(matrix, syndrome, syndrome),
        )
        if result.status not in (0, 1, 2):
            raise RuntimeError(f"MLE optimizer failed: {result.message}")
        return (None if result.x is None else result.x[:n],
                {"optimal": result.status == 0,
                 "status": ("optimal" if result.status == 0 else
                            "infeasible" if result.status == 2 else "interrupted")})

    def _copt(self, syndrome):
        try:
            import coptpy as cp
            from coptpy import COPT
        except ImportError as exc:
            raise ImportError("Install coptpy and configure a COPT license") from exc
        m, n = self.h.shape
        env = model = None
        try:
            env = cp.Envr()
            model = env.createModel("dem_mle")
            model.setParam("Logging", 0)
            matrix = hstack((self.h, -2 * eye(m, format="csc")), format="csc")
            model.loadMatrix(
                c=np.concatenate((self.weights, np.zeros(m))), A=matrix,
                lhs=syndrome.astype(float), rhs=syndrome.astype(float),
                lb=np.zeros(n + m),
                ub=np.concatenate((np.ones(n), np.floor(self.row_weights / 2))),
                vtype=[COPT.BINARY] * n + [COPT.INTEGER] * m,
            )
            model.ObjSense = COPT.MINIMIZE
            model.solve()
            status = model.MipStatus
            if status in (COPT.NUMERICAL, COPT.UNBOUNDED, COPT.INF_OR_UNB):
                raise RuntimeError(f"MLE optimizer failed with status {status}")
            variables = model.getVars()
            values = (np.asarray(model.getInfo(
                COPT.Info.Value, [variables[j] for j in range(n)]))
                if model.HasMipSol else None)
            return values, {
                "optimal": status == COPT.OPTIMAL,
                "status": ("optimal" if status == COPT.OPTIMAL else
                           "infeasible" if status == COPT.INFEASIBLE else "interrupted"),
            }
        finally:
            if model is not None:
                model.dispose()
            if env is not None:
                env.dispose()

    def decode(self, syndrome):
        """Return a syndrome-checked correction and optimizer diagnostics.

        ``accepted`` denotes a feasible correction, not logical correctness.
        ``diagnostics['optimal']`` denotes solver optimal termination under its
        numerical tolerances. Backend exceptions are never counted as success.
        """
        syndrome = _binary(syndrome, self.h.shape[0], "Syndrome")
        started = time.perf_counter()
        if not syndrome.any():
            correction = np.zeros(self.h.shape[1], dtype=np.uint8)
            details = {"status": "optimal", "optimal": True}
        else:
            raw, details = getattr(self, "_" + self.backend)(syndrome)
            correction = None
            if raw is not None:
                rounded = np.rint(np.asarray(raw))
                if not np.all(np.isfinite(raw)) or not np.allclose(
                        raw, rounded, rtol=0, atol=1e-5):
                    raise RuntimeError("Optimizer returned a nonintegral correction")
                correction = _binary(rounded, self.h.shape[1], "Optimizer correction")
                if not _valid(self.h, correction, syndrome):
                    raise RuntimeError("Optimizer correction failed the exact parity check")
        accepted = correction is not None
        details.update(decoder="mle", backend=self.backend,
                       objective=float(self.weights @ correction) if accepted else None,
                       seconds=time.perf_counter() - started)
        return {"correction": correction, "accepted": accepted,
                "diagnostics": details}
