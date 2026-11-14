"""Physical CXSWAP memory circuits for the Armillary candidate catalog."""
from __future__ import annotations

import hashlib

import numpy as np
import stim

from .algebra import indices, legacy_logicals


def build_circuit(record, basis: str, rounds: int = 4, p: float = 0.001):
    """Construct an X- or Z-memory experiment in fixed physical-site coordinates.

    The first and last syndrome sweeps are ideal. Between them, ``rounds``
    noisy sweeps alternate forward and reverse transport. The returned audit
    includes a small noiseless consistency check and a canonical circuit hash.
    """
    if basis not in ("x", "z"):
        raise ValueError("basis must be 'x' or 'z'")
    if isinstance(rounds, bool) or not isinstance(rounds, int) or rounds < 0:
        raise ValueError("rounds must be a nonnegative integer")
    if not np.isfinite(p) or not 0 <= p < 0.5:
        raise ValueError("p must be a finite probability in [0, 0.5)")
    n, period, q = record["n"], record["p"], record["q"]
    P = q * period
    if n != 12 * P or n % 4:
        raise ValueError("unsupported q3/q4 block dimensions")
    hx, hz = record["_hx"], record["_hz"]

    logical_constraints = hz if basis == "x" else hx
    logical_stabilizers = hx if basis == "x" else hz
    checks = hx if basis == "x" else hz
    logicals = legacy_logicals(logical_constraints, logical_stabilizers, n)
    if len(logicals) != record["k"]:
        raise ValueError(f"logical basis has {len(logicals)} rows, expected {record['k']}")

    # Track role labels through the physical endpoint swaps, then emit the
    # Stim circuit in fixed physical-site coordinates.
    label_site = {}
    site_label = {}
    for typ in (0, 1):
        for group in range(12):
            for z in range(period):
                for lane in range(q):
                    label = typ * n + group * P + z * q + lane
                    site = (typ, (-group) % 12, z, lane)
                    label_site[label] = site
                    site_label[site] = label
    initial = label_site.copy()
    site_id = {site: i for i, site in enumerate(sorted(site_label))}
    top, phases = record["top"], record["phases"]
    mode = record["mode"]
    x_groups = set(record["skeleton"]["xrows"])
    z_groups = set(record["skeleton"]["zrows"])
    routes = []
    for t in range(12):
        pairs = []
        for group in range(12):
            for member in range(P):
                anc = n + group * P + member
                typ, c, z, lane = label_site[anc]
                if t % 2 == 0:
                    if isinstance(top, dict):
                        lane2 = top["permutations"][f"{t},{group}"][lane]
                    else:
                        signs, shifts = top
                        lane2 = ((-1) ** signs[t * 12 + group] * (lane - shifts[t * 12 + group])) % q
                    dst = (0, c, z, lane2)
                else:
                    j = ((t // 2) * 12 + group) * (q if mode == "channel" else 1) + (lane if mode == "channel" else 0)
                    dst = (1, (c - 1) % 12, (z + phases[j]) % period, lane)
                data = site_label[dst]
                if data >= n:
                    raise ValueError("CXSWAP route did not meet a data-role site")
                kind = "x" if group in x_groups else "z" if group in z_groups else None
                pairs.append((anc, data, label_site[anc], label_site[data], kind))
        touched = [v for edge in pairs for v in edge[:2]]
        if len(set(touched)) != 2 * n:
            raise ValueError(f"physical layer {t} is not a full matching")
        for a, b, _, _, _ in pairs:
            sa, sb = label_site[a], label_site[b]
            label_site[a], label_site[b] = sb, sa
            site_label[sa], site_label[sb] = b, a
        routes.append(pairs)
    final = label_site.copy()

    def locs(mapping, labels):
        return [site_id[mapping[i]] for i in labels]

    x_checks = [n + g * P + i for g in record["skeleton"]["xrows"] for i in range(P)]
    z_checks = [n + g * P + i for g in record["skeleton"]["zrows"] for i in range(P)]
    cx_check_labels = x_checks
    half = n // 4
    if len(x_checks) != half or len(z_checks) != half:
        raise ValueError("unexpected q3/q4 check count")

    circuit = stim.Circuit()
    initial_data = locs(initial, range(n))
    circuit.append("R", initial_data)
    if basis == "x":
        circuit.append("H", initial_data)
    for r in range(rounds + 2):
        forward = (r % 2) == 0
        before, after = (initial, final) if forward else (final, initial)
        noisy = 0 < r <= rounds
        anc = locs(before, range(n, 2 * n))
        data = locs(before, range(n))
        circuit.append("R", anc)
        if noisy:
            # Preparation noise on auxiliary transport roles, before the
            # ideal basis rotation for X checks.
            circuit.append("DEPOLARIZE1", anc, p)
        check_h = locs(before, cx_check_labels)
        circuit.append("H", check_h)
        circuit.append("TICK")
        for layer in routes if forward else reversed(routes):
            cx, swaps = [], []
            for a, b, src_site, dst_site, kind in layer:
                su, sv = site_id[src_site], site_id[dst_site]
                u, v = (su, sv) if forward else (sv, su)
                swaps.extend((u, v))
                if kind == "x":
                    cx.extend((u, v))
                elif kind == "z":
                    cx.extend((v, u))
            if cx:
                circuit.append("CX", cx)
            circuit.append("SWAP", swaps)
            if noisy:
                circuit.append("DEPOLARIZE2", swaps, p)
            circuit.append("TICK")
        check_h = locs(after, cx_check_labels)
        circuit.append("H", check_h)
        xloc, zloc = locs(after, x_checks), locs(after, z_checks)
        if noisy:
            # Measurement noise after the ideal basis rotation. Data have only
            # this one idle interval, representing measurement plus preparation.
            circuit.append("DEPOLARIZE1", xloc + zloc, p)
            circuit.append("DEPOLARIZE1", locs(after, range(n)), p)
        circuit.append("M", xloc + zloc)
        offset = 0 if basis == "x" else half
        for j in range(half):
            refs = [stim.target_rec(-2 * half + offset + j)]
            if r:
                refs.append(stim.target_rec(-4 * half + offset + j))
            circuit.append("DETECTOR", refs)
        circuit.append("TICK")
    last_after = initial if (rounds + 1) % 2 == 1 else final
    data = locs(last_after, range(n))
    if basis == "x":
        circuit.append("H", data)
    circuit.append("M", data)
    offset = 0 if basis == "x" else half
    for j, row in enumerate(checks):
        circuit.append("DETECTOR", [stim.target_rec(-n + v) for v in indices(row)] + [stim.target_rec(-n - 2 * half + offset + j)])
    for j, row in enumerate(logicals):
        circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-n + v) for v in indices(row)], j)
    noiseless = circuit.without_noise().compile_detector_sampler(seed=261102)
    dets, obs = noiseless.sample(32, separate_observables=True)
    noiseless_pass = not (np.any(dets) or np.any(obs))
    return circuit, {
        "matrix_sha256": record["matrix_sha256"], "n": n, "k": record["k"],
        "basis": basis, "rounds": rounds, "p": p, "detectors": circuit.num_detectors,
        "observables": circuit.num_observables,
        "noiseless_pass": bool(noiseless_pass),
        "noiseless_detector_events": int(dets.sum()),
        "noiseless_logical_events": int(obs.sum()),
        "circuit_sha256": hashlib.sha256(str(circuit).encode()).hexdigest(),
        "candidate_id": record["candidate_id"],
    }
