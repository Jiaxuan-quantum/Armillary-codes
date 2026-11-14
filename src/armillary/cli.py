"""Command-line entry points for the reference workflow."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(prog="armillary", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    catalog = commands.add_parser("catalog", help="List the reproducible candidate catalog")
    catalog.add_argument("--distance", type=int)
    verify = commands.add_parser("verify", help="Reconstruct matrices and check canonical hashes")
    verify.add_argument("candidate", nargs="?", default="all")
    model = commands.add_parser("model", help="Build and save a memory circuit and detector model")
    model.add_argument("candidate")
    model.add_argument("--basis", choices=("x", "z"), required=True)
    model.add_argument("--p", type=float, default=.001)
    model.add_argument("--rounds", type=int, default=4)
    model.add_argument("--out", type=Path, required=True)
    sample = commands.add_parser("sample", help="Sample a saved circuit and run short-to-long Relay")
    sample.add_argument("model", type=Path)
    sample.add_argument("--shots", type=int, default=32)
    sample.add_argument("--seed", type=int, required=True)
    sample.add_argument("--batch-size", type=int, default=32)
    sample.add_argument("--out", type=Path, required=True)
    decode = commands.add_parser("decode", help="Decode saved Relay handoffs with MLE")
    decode.add_argument("model", type=Path)
    decode.add_argument("front_run", type=Path)
    decode.add_argument("--backend", choices=("scipy", "copt"), default="scipy")
    decode.add_argument("--out", type=Path, required=True)
    commands.add_parser("benchmark", help="Recompute the released aggregate benchmark estimates")
    args = parser.parse_args(argv)
    if args.command == "catalog":
        from .candidates import list_candidates
        rows = list_candidates()
        if args.distance is not None:
            rows = [r for r in rows if r["exact_distance"] == args.distance]
        print(json.dumps(rows, indent=2))
    elif args.command == "verify":
        from .candidates import list_candidates, load_candidate
        names = [r["candidate_id"] for r in list_candidates()] if args.candidate == "all" else [args.candidate]
        for name in names:
            record = load_candidate(name)
            print(f"PASS {record['candidate_id']} n={record['n']} k={record['k']} "
                  f"sha256={record['matrix_sha256']}", flush=True)
        print(f"Verified {len(names)} candidate(s). Distance searches were not rerun.")
    elif args.command == "model":
        from .models import build_model
        result = build_model(args.candidate, args.basis, rounds=args.rounds, p=args.p)
        result.save(args.out)
        print(json.dumps(result.metadata, indent=2))
    elif args.command == "sample":
        from .experiments import sample
        sample(args.model, args.shots, args.seed, args.out, args.batch_size)
    elif args.command == "decode":
        from .experiments import decode_handoffs
        decode_handoffs(args.model, args.front_run, args.out, backend=args.backend)
    elif args.command == "benchmark":
        from importlib.resources import files
        from .statistics import failure_statistics
        data = json.loads(files("armillary").joinpath("data/benchmark.json").read_text(encoding="utf-8"))
        for row in data["pooled"]:
            statistics = failure_statistics(row["W1"] + row["W2"], row["N"], row["k"], data["rounds"])
            print(json.dumps(dict(candidate=row["candidate"], **statistics)))


if __name__ == "__main__":
    main()
