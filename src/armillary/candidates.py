"""Load and reconstruct the finite-displacement routing-code catalog.

The matrices are recovered by replaying the physical swaps and CNOT actions,
not by loading precomputed parity-check arrays. Candidate distance values are
reported metadata; reconstruction does not perform a minimum-distance search.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from importlib.resources import files
import hashlib
import json

from .algebra import echelon


def list_candidates() -> list[dict]:
    """Return the compact catalog, without reconstructing any matrices."""
    source = files("armillary").joinpath("data", "catalog.json")
    return json.loads(source.read_text(encoding="utf-8"))["candidates"]


def candidate_definition(identifier: str) -> dict:
    """Read a definition by full candidate ID or one of its short aliases."""
    identifier = str(identifier).lower()
    matches = [row for row in list_candidates()
               if identifier in [row["candidate_id"].lower(),
                                 *[s.lower() for s in row.get("aliases", [])]]]
    if len(matches) != 1:
        raise ValueError(f"Unknown or ambiguous candidate: {identifier}")
    name = matches[0]["candidate_id"] + ".json"
    source = files("armillary").joinpath("data", "candidates", name)
    encoded = source.read_bytes()
    if hashlib.sha256(encoded).hexdigest() != matches[0]["definition_sha256"]:
        raise ValueError("Candidate definition checksum disagrees with catalog")
    return json.loads(encoded)


def _permutation(record, step: int, group: int) -> list[int]:
    q = record["q"]
    top = record["top"]
    if isinstance(top, dict):
        return top["permutations"][f"{step},{group}"]
    signs, shifts = top
    offset = step * 12 + group
    return [((-1) ** signs[offset] * (j - shifts[offset])) % q for j in range(q)]


def matrix_digest(n: int, hx: list[int], hz: list[int]) -> str:
    """Hash the ordered matrices in their canonical hexadecimal-row format."""
    encoded = json.dumps([n, list(map(hex, hx)), list(map(hex, hz))],
                         separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def reconstruct(record: dict) -> dict:
    """Replay one forward and one reverse sweep and return the CSS matrices.

    All twelve ancilla groups participate in transport, including groups whose
    stabilizer checks were removed. Only retained X/Z groups apply CNOTs.
    Bit propagation checks that each sweep extracts its intended parity checks
    without residual ancilla coupling. The reverse sweep must restore all slots
    and reproduce exactly the same ordered check supports.
    """
    period, q = record["p"], record["q"]
    if q not in (3, 4) or not isinstance(period, int) or period <= 0:
        raise ValueError("Expected q in {3, 4} and a positive integer period")
    channel = record["mode"] == "channel"
    if record["mode"] not in ("common", "channel"):
        raise ValueError("Unknown displacement mode")
    phases = record["phases"]
    if len(phases) != 72 * (q if channel else 1):
        raise ValueError("Incorrect number of displacement parameters")
    for step in range(0, 12, 2):
        for group in range(12):
            if sorted(_permutation(record, step, group)) != list(range(q)):
                raise ValueError("Invalid lane permutation")
    block = q * period
    n = 12 * block
    masks = [record["skeleton"]["xrows"], record["skeleton"]["zrows"]]
    xs, zs = map(set, masks)
    if (len(xs) != 3 or len(zs) != 3 or xs & zs
            or not (xs | zs) <= set(range(12))):
        raise ValueError("Expected three distinct X and three distinct Z groups")

    def site(rail, column, z, lane):
        return (rail * 12 + column) * block + q * z + lane

    # Ancilla labels occupy [0,n); data labels occupy [n,2n).
    at, position = [0] * (2 * n), [0] * (2 * n)
    for group in range(12):
        for z in range(period):
            for lane in range(q):
                ancilla = group * block + q * z + lane
                data = n + ancilla
                ancilla_site = site(1, (-group) % 12, z, lane)
                data_site = site(0, (-group) % 12, z, lane)
                at[ancilla_site], position[ancilla] = ancilla, ancilla_site
                at[data_site], position[data] = data, data_site
    initial = at[:]
    steps = []
    for step in range(12):
        pairs, covered = [], set()
        for group in range(12):
            permutation = _permutation(record, step, group) if step % 2 == 0 else None
            for member in range(block):
                label = group * block + member
                source = position[label]
                site_group, member_at_site = divmod(source, block)
                rail, column = divmod(site_group, 12)
                z, lane = divmod(member_at_site, q)
                if step % 2 == 0:
                    if rail != 1:
                        raise ValueError("Unexpected ancilla rail in an even layer")
                    destination = site(0, column, z, permutation[lane])
                else:
                    if rail != 0:
                        raise ValueError("Unexpected ancilla rail in an odd layer")
                    index = (step // 2 * 12 + group) * (q if channel else 1)
                    index += lane if channel else 0
                    destination = site(1, (column - 1) % 12,
                                       (z + phases[index]) % period, lane)
                if source in covered or destination in covered or at[destination] < n:
                    raise ValueError("A transport layer is not an ancilla-data matching")
                covered.update((source, destination))
                pairs.append((source, destination))
        if len(covered) != 2 * n:
            raise ValueError("A transport layer does not cover every physical site")
        for u, v in pairs:
            a, b = at[u], at[v]
            at[u], at[v], position[a], position[b] = b, a, v, u
        steps.append(pairs)

    at = initial[:]
    sweeps = []
    for sequence in (steps, list(reversed(steps))):
        bits = [1 << label for label in at]
        visited = defaultdict(list)
        for pairs in sequence:
            for u, v in pairs:
                ancilla_site, data_site = (u, v) if at[u] < n else (v, u)
                ancilla, data = at[ancilla_site], at[data_site] - n
                if not (0 <= ancilla < n and 0 <= data < n):
                    raise ValueError("Unexpected role labels during a sweep")
                visited[ancilla].append(data)
                if ancilla // block in xs:
                    bits[data_site] ^= bits[ancilla_site]
                elif ancilla // block in zs:
                    bits[ancilla_site] ^= bits[data_site]
                at[u], at[v] = at[v], at[u]
                bits[u], bits[v] = bits[v], bits[u]
        incidence = defaultdict(int)
        for ancilla, support in visited.items():
            if (len(support) != 12 or len(set(support)) != 12
                    or sorted(data // block for data in support) != list(range(12))):
                raise ValueError("An ancilla must meet one data qubit in each group")
            if ancilla // block in xs:
                for data in support:
                    incidence[data] ^= 1 << ancilla
        for slot, label in enumerate(at):
            expected = 1 << label
            if label >= n:
                expected ^= incidence[label - n]
            elif label // block in zs:
                expected ^= sum(1 << (n + data) for data in visited[label])
            if bits[slot] != expected:
                raise ValueError("A sweep contains an unwanted residual CNOT action")
        sweeps.append([[sorted(visited[group * block + j])
                        for group in mask for j in range(block)] for mask in masks])
    if at != initial or sweeps[0] != sweeps[1]:
        raise ValueError("Forward/reverse transport or parity checks disagree")
    hx, hz = [[sum(1 << data for data in row) for row in rows] for rows in sweeps[0]]
    # Two common columns form a Tanner-graph 4-cycle.
    girth_ge6 = all(max(Counter((row[i], row[j]) for row in rows
                               for i in range(12) for j in range(i + 1, 12)).values()) == 1
                   for rows in sweeps[0])
    rank_x, rank_z = len(echelon(hx)), len(echelon(hz))
    return {"n": n, "k": n - rank_x - rank_z, "hx": hx, "hz": hz,
            "rank_x": rank_x, "rank_z": rank_z,
            "matrix_sha256": matrix_digest(n, hx, hz),
            "girth_ge6": girth_ge6, "slots_restored": True,
            "forward_reverse_agree": True}


def load_candidate(identifier: str) -> dict:
    """Load a catalog definition, reconstruct it, and check its identity."""
    record = candidate_definition(identifier)
    rebuilt = reconstruct(record)
    for key in ("n", "k", "matrix_sha256"):
        if record[key] != rebuilt[key]:
            raise ValueError(f"Candidate reconstruction disagrees with recorded {key}")
    record["_hx"], record["_hz"] = rebuilt["hx"], rebuilt["hz"]
    record["_verification"] = {key: value for key, value in rebuilt.items()
                               if key not in ("hx", "hz")}
    return record
