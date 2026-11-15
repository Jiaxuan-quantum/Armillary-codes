"""Reconstruct and verify the three published continuous hardware layouts.

This entry point is independent of the Armillary installation and can be run
from any working directory. It requires Python 3.11+ and SciPy. It does not
perform a new path search, distance search, or quantum-error-correction run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import platform
import sys
import time


HERE = Path(__file__).resolve().parent
CODE_NAMES = ("d12", "d14", "d16")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8", newline="\n")


def check_manifest():
    """Check the release inputs before importing verification modules."""
    manifest = read_json(HERE / "MANIFEST.json")
    require(manifest.get("schema") == "armillary-layout-manifest-v1",
            "Unsupported manifest schema")
    entries = manifest.get("sha256")
    require(isinstance(entries, dict) and entries, "Empty manifest")
    expected = {"reproduce.py", "reference_results.json", "PROVENANCE.json",
                "requirements.txt", "README.md"}
    for folder in ("modules", "data"):
        expected.update(p.relative_to(HERE).as_posix()
                        for p in (HERE / folder).rglob("*")
                        if p.is_file() and p.suffix in (".py", ".json"))
    require(set(entries) == expected, "Manifest does not cover the release input set")
    for name, digest in entries.items():
        relative = PurePosixPath(name)
        require(not relative.is_absolute() and ".." not in relative.parts
                and "\\" not in name and ":" not in name,
                f"Unsafe manifest path: {name}")
        path = (HERE / relative).resolve()
        require(path.is_relative_to(HERE), f"Manifest path leaves the package: {name}")
        require(path.is_file(), f"Missing release input: {name}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        require(actual == digest, f"Checksum mismatch: {name}")
    return len(entries), hashlib.sha256((HERE / "MANIFEST.json").read_bytes()).hexdigest()


def close(actual, expected, label, abs_tol=1e-8):
    require(math.isclose(actual, expected, abs_tol=abs_tol, rel_tol=1e-11),
            f"Mismatch for {label}: {actual} != {expected}")


def verify_one(tag, output, reference, api):
    started = time.perf_counter()
    source = HERE / "data" / tag
    geometry = read_json(source / "geometry.json")
    hardware = read_json(source / "hardware.json")
    candidate = read_json(source / "candidate.json")
    stack = read_json(source / "stack.json")
    require(geometry["code"] == reference["code"], f"Wrong geometry identity for {tag}")
    require(candidate["candidate_id"] == reference["candidate_id"],
            f"Wrong candidate identity for {tag}")
    require([candidate["n"], candidate["k"], candidate["exact_distance"]]
            == reference["code"], f"Wrong candidate parameters for {tag}")

    model = api["build"](geometry, hardware)
    frozen = read_json(source / "routes.json")
    # The d16 source uses an earlier descriptive schema string. Its numerical
    # representation is identical. Compare every numerical/discrete field.
    ignored = {"schema", "interpretation"}
    require(set(model) == set(frozen), f"Route fields differ for {tag}")
    for key in set(model) - ignored - {"orbits"}:
        require(model[key] == frozen[key], f"Route metadata differs: {tag}/{key}")
    require(len(model["orbits"]) == len(frozen["orbits"]), "Route count mismatch")
    for actual, expected in zip(model["orbits"], frozen["orbits"]):
        require(set(actual) == set(expected), "Route record fields differ")
        for key in ("key", "face", "source", "target", "points_ur", "endpoints_ur"):
            require(actual[key] == expected[key], f"Route differs: {tag}/{actual['key']}/{key}")
        require(len(actual["segment_arc_lengths"]) == len(expected["segment_arc_lengths"]),
                "Segment count differs")
        for a, b in zip(actual["segment_arc_lengths"], expected["segment_arc_lengths"]):
            close(a, b, "segment arc length", 1e-9)
        close(actual["length"], expected["length"], "connection arc length", 1e-9)
        close(actual["input_length"], expected["input_length"], "stored connection length", 1e-9)
    max_length_error = max(abs(e["length"] - e["input_length"]) for e in model["orbits"])
    require(max_length_error < 1e-8, "Integrated and input arc lengths differ")

    code = api["literal_code"](candidate, geometry["q"])
    require(code["matrix_sha256"] == candidate["matrix_sha256"] == reference["matrix_sha256"],
            f"Matrix identity mismatch for {tag}")
    require([code["n"], code["k"]] == reference["code"][:2], "Code dimension mismatch")
    require(code["girth_ge6"] and code["native_both_rounds_literal"]
            and code["slots_restored"], "Code/circuit replay failed")
    witness_weight = api["verify_logical"](code, read_json(source / "logical_witness.json"))
    require(witness_weight == reference["code"][2], "Logical witness has wrong weight")

    graph = read_json(source / "conflict_graph.json")
    colors = [e["layer"] for e in hardware["edge_orbits"]]
    require(graph["keys"] == [e["key"] for e in hardware["edge_orbits"]], "Graph key mismatch")
    require(not graph["self_conflicts"], "Self-conflicting route")
    require(all(colors[i] != colors[j] for i, j in graph["edges"]),
            "The selected faces violate a conflict edge")

    checks = {
        "continuous": api["audit"](source / "geometry.json", source / "hardware.json"),
        "route_export": api["compare_model"](source / "geometry.json", source / "hardware.json", model),
        "circuit": api["replay"](candidate, model),
    }
    for name, check in checks.items():
        require(check["status"] == "PASS", f"{tag}/{name} failed")
        write_json(output / f"{name}.json", check)

    metrics = api["stats"](model, stack["tier_faces"])
    for name, expected in reference["metrics"].items():
        if isinstance(expected, (int, float)):
            close(metrics[name], expected, f"{tag}/{name}")
    require(metrics["max_face_endpoint_degree"] <= 4, "Per-face endpoint degree exceeds four")
    for face in set(colors):
        edges = [e for e in model["orbits"] if e["face"] == face]
        require(hardware["layer_counts"][str(face)] == len(edges) * model["p"],
                "Face population mismatch")
        close(hardware["layer_length_sums"][str(face)],
              sum(e["length"] for e in edges) * model["p"], "Face length sum", 1e-7)
    exhaustive = api["optimize"](hardware, source=f"data/{tag}/hardware.json")
    close(exhaustive["minimum_cost"]["C_hw"], metrics["C_hw"], "Fixed-face optimum", 1e-11)
    require(exhaustive["number_of_partitions"] == stack["number_of_partitions"]
            == reference["stack_partitions_checked"], "Stacking enumeration count mismatch")

    wire = checks["continuous"]["minimum_wire_lower_bound"]
    device = checks["continuous"]["minimum_device_region_lower_bound"]
    require(wire > 0.1 and device > 0.05, "A strict geometric requirement failed")
    close(wire, reference["minimum_wire_clearance_lower_bound"], "Wire separation bound", 1e-10)
    close(device, reference["minimum_nonendpoint_device_region_clearance_lower_bound"],
          "Device-region clearance bound", 1e-10)
    close(checks["circuit"]["minimum_device_center_distance"], 1.0, "Minimum qubit pitch", 1e-10)
    result = dict(
        status="PASS", code=reference["code"], candidate_id=candidate["candidate_id"],
        matrix_sha256=code["matrix_sha256"], metrics=metrics,
        minimum_wire_clearance_lower_bound=wire,
        minimum_nonendpoint_device_region_clearance_lower_bound=device,
        minimum_qubit_center_distance=checks["circuit"]["minimum_device_center_distance"],
        maximum_integrated_length_discrepancy=max_length_error,
        route_orbits=len(model["orbits"]), connections=model["total_connections"],
        positions=model["total_devices"], graph_constraints_checked=len(graph["edges"]),
        stack_partitions_checked=exhaustive["number_of_partitions"],
        matrix_reconstruction={k: v for k, v in code.items() if k not in ("hx", "hz")},
        logical_witness_weight=witness_weight, distance_scope=reference["distance_scope"],
        seconds=time.perf_counter() - started,
    )
    write_json(output / "routes.json", model)
    write_json(output / "METRICS.json", metrics)
    write_json(output / "STACKING.json", {
        "status": exhaustive["status"],
        "number_of_partitions": exhaustive["number_of_partitions"],
        "minimum_cost": exhaustive["minimum_cost"], "scope": exhaustive["scope"],
    })
    write_json(output / "VERIFY.json", result)
    print(json.dumps({"code": tag, "status": "PASS", "tiers": metrics["tiers"],
                      "C_hw": metrics["C_hw"], "seconds": result["seconds"]}), flush=True)
    return result


def main():
    if not __debug__:
        raise SystemExit("Verification requires assertions. Run without Python -O or PYTHONOPTIMIZE.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code", choices=(*CODE_NAMES, "all"), default="all")
    parser.add_argument("--out", type=Path, required=True,
                        help="New or empty output directory outside hardware_layouts")
    args = parser.parse_args()
    started = time.perf_counter()
    count, manifest_sha = check_manifest()
    output = args.out.resolve()
    require(not output.is_relative_to(HERE), "Output must be outside the immutable input package")
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())),
            "Output directory must be new or empty")
    reference = read_json(HERE / "reference_results.json")
    sys.path.insert(0, str(HERE / "modules"))
    from route_model import build, stats
    from native import literal_code, verify_logical
    from audit_replay import replay
    from continuous_audit import audit, compare_model
    from optimize_stack import optimize
    import scipy

    api = dict(build=build, stats=stats, literal_code=literal_code,
               verify_logical=verify_logical, replay=replay, audit=audit,
               compare_model=compare_model, optimize=optimize)
    selected = CODE_NAMES if args.code == "all" else (args.code,)
    output.mkdir(parents=True, exist_ok=True)
    summary = dict(
        status="RUNNING", codes={}, manifest_files_checked=count,
        manifest_sha256=manifest_sha, python=platform.python_version(), scipy=scipy.__version__,
        metric_scope=reference["metric_scope"],
        verification_scope="Complete fixed-edge circuit replay, continuous geometry, and exhaustive stacking for fixed faces. Paths and face assignments are not globally reoptimized.",
        limitations="Double-precision geometric bounds, not a formal interval proof or manufacturing acceptance. No shielding, vertical-via footprint, control/readout routing, electromagnetic simulation, new distance search, or QEC sampling.",
    )
    write_json(output / "SUMMARY.json", summary)
    try:
        for tag in selected:
            folder = output / tag
            folder.mkdir()
            summary["codes"][tag] = verify_one(tag, folder, reference["codes"][tag], api)
            write_json(output / "SUMMARY.json", summary)
    except Exception as exc:
        summary.update(status="FAILED", error=f"{type(exc).__name__}: {exc}",
                       seconds=time.perf_counter() - started)
        write_json(output / "SUMMARY.json", summary)
        raise
    summary.update(status="PASS", seconds=time.perf_counter() - started)
    write_json(output / "SUMMARY.json", summary)
    print(json.dumps({"status": "PASS", "codes": list(summary["codes"]),
                      "manifest_files_checked": count, "seconds": summary["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
