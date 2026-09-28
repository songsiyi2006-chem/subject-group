"""Reproduce the lightweight HVI differential check used in report review.

This compares two separate implementations in metrics_control.py: vectorized
single-candidate area subtraction and explicit before/after rectangle unions.
It is not an external-library validation or a proof for every possible input.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform

import numpy as np

from metrics_control import hypervolume_2d, hypervolume_improvement


BASE = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    seed, batches, observed_count, candidates_per_batch = 718, 20, 10, 20
    rng = np.random.default_rng(seed)
    reference = np.array([-0.5, -0.75])
    rows = []
    for batch in range(batches):
        observed = rng.normal(size=(observed_count, 2))
        candidates = rng.normal(size=(candidates_per_batch, 2))
        baseline = hypervolume_2d(observed, reference)
        explicit = np.array([
            hypervolume_2d(np.vstack([observed, candidate]), reference) - baseline
            for candidate in candidates
        ])
        vectorized = hypervolume_improvement(candidates, observed, reference)
        rows.append({
            "batch_index": batch,
            "observed_vectors": observed.tolist(),
            "candidate_draws": candidates.tolist(),
            "explicit_union_difference": explicit.tolist(),
            "vectorized_single_candidate_HVI": vectorized.tolist(),
            "maximum_absolute_difference": float(np.max(np.abs(explicit - vectorized))),
        })
    maximum = max(row["maximum_absolute_difference"] for row in rows)
    tolerance = 1e-12
    result = {
        "schema_version": "electratwin.cross_review.v1",
        "evidence_role": "executed_numerical_differential_check_not_experimental_validation",
        "executed_utc": datetime.now(timezone.utc).isoformat(),
        "random_seed": seed,
        "generator": "numpy.random.default_rng; PCG64",
        "distribution": "Independent standard normal N(0,1) coordinates for each observed and candidate vector; a fresh observed set is drawn before the candidate set in each batch.",
        "batches": batches,
        "observed_vectors_per_batch": observed_count,
        "candidate_draws_per_batch": candidates_per_batch,
        "total_candidate_comparisons": batches * candidates_per_batch,
        "reference_point": reference.tolist(),
        "comparison": "hypervolume_improvement(draw, observed, reference) versus hypervolume_2d(observed union draw, reference) minus hypervolume_2d(observed, reference)",
        "maximum_absolute_difference": maximum,
        "absolute_tolerance": tolerance,
        "passed": bool(maximum < tolerance),
        "scope": [
            "Reproduces the report-review calculation; not an additional unit test counted in the 36-test suite.",
            "Differential comparison of two distinct geometry implementations within the same module; not independent external-library or formal-method validation.",
            "Tests raw HVI geometry, including vectors below the reference; does not test GP calibration, acquisition sampling quality, physical-domain masking, chemistry or instrument operation.",
            "No PDE campaign, network access, physical instrument connection or experimental evaluation is performed.",
        ],
        "source_sha256": {
            "electratwin/scripts/cross_review.py": sha256(Path(__file__)),
            "electratwin/scripts/metrics_control.py": sha256(BASE / "scripts" / "metrics_control.py"),
        },
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "scipy": version("scipy"), "scikit_learn": version("scikit-learn")},
        "batches_data": rows,
    }
    output = BASE / "results" / "cross_review.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"total_candidate_comparisons": result["total_candidate_comparisons"],
                      "maximum_absolute_difference": maximum, "passed": result["passed"]}))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
