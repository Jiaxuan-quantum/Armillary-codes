"""Release boundaries for the standalone geometric reproduction inputs."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
LAYOUTS = ROOT / "hardware_layouts"


def loader():
    spec = importlib.util.spec_from_file_location("layout_reproduction", LAYOUTS / "reproduce.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_layout_candidates_match_public_catalog():
    reference = json.loads((LAYOUTS / "reference_results.json").read_text(encoding="utf-8"))
    for tag, entry in reference["codes"].items():
        candidate = LAYOUTS / "data" / tag / "candidate.json"
        published = ROOT / "src/armillary/data/candidates" / (entry["candidate_id"] + ".json")
        assert candidate.read_bytes() == published.read_bytes()
        record = json.loads(candidate.read_text(encoding="utf-8"))
        assert record["matrix_sha256"] == entry["matrix_sha256"]
        assert [record["n"], record["k"], record["exact_distance"]] == entry["code"]


def test_layout_manifest_rejects_modified_input(tmp_path):
    module = loader()
    count, _ = module.check_manifest()
    assert count > 0
    destination = tmp_path / "hardware_layouts"
    shutil.copytree(LAYOUTS, destination, ignore=shutil.ignore_patterns("__pycache__"))
    module.HERE = destination
    target = destination / "data/d16/geometry.json"
    target.write_bytes(target.read_bytes() + b" ")
    with pytest.raises(ValueError, match="Checksum mismatch: data/d16/geometry.json"):
        module.check_manifest()


def test_layout_reproduction_rejects_optimized_python(tmp_path):
    result = subprocess.run(
        [sys.executable, "-O", str(LAYOUTS / "reproduce.py"), "--out", str(tmp_path / "output")],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "requires assertions" in result.stderr
    assert not (tmp_path / "output").exists()


def test_layout_reproduction_rejects_output_inside_inputs():
    result = subprocess.run(
        [sys.executable, str(LAYOUTS / "reproduce.py"), "--out", str(LAYOUTS / "forbidden_output")],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "outside the immutable input package" in result.stderr
    assert not (LAYOUTS / "forbidden_output").exists()
