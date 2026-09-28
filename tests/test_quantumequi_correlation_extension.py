"""Fast invariants and independent saved-output audits; no Psi4 imports or jobs."""
import csv
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "quantumequi/scripts/correlation_extension.py"
spec = importlib.util.spec_from_file_location("correlation_extension", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
OUT = REPO / "quantumequi/results/extensions/correlation"


class AnalyticChecks(unittest.TestCase):
    def test_paired_singlet_and_doublet_spin(self):
        self.assertAlmostEqual(m.spin_squared_from_occupied(np.eye(2)[:, :1], np.eye(2)[:, :1], np.eye(2)), 0)
        self.assertAlmostEqual(m.spin_squared_from_occupied(np.eye(2)[:, :1], np.empty((2, 0)), np.eye(2)), .75)

    def test_separated_spin_and_common_rotation(self):
        a, b = np.eye(2)[:, :1], np.eye(2)[:, 1:]
        theta = .719
        q = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        self.assertAlmostEqual(m.spin_squared_from_occupied(a, b, np.eye(2)), 1)
        self.assertAlmostEqual(m.spin_squared_from_occupied(q @ a, q @ b, np.eye(2)), 1)

    def test_invalid_orbitals_rejected(self):
        with self.assertRaises(ValueError):
            m.spin_squared_from_occupied(np.ones((2, 1)), np.ones((3, 1)), np.eye(2))

    def test_five_point_quartic_derivatives(self):
        h = .025
        x = np.arange(-2, 3) * h
        values = 3 + 2 * x + 7 * x**2 + 11 * x**3 + 13 * x**4
        found = m.five_point_derivatives(values, h)
        self.assertAlmostEqual(found["gradient_5point_Hartree_A"], 2, places=11)
        self.assertAlmostEqual(found["curvature_Hartree_A2"], 14, places=10)
        with self.assertRaises(ValueError):
            m.five_point_derivatives(values, 0)

    def test_branch_requires_energy_and_spin(self):
        self.assertEqual(m.classify_uhf(-1, -.8, .9), "observed_lower_spin_broken_branch")
        self.assertEqual(m.classify_uhf(-1, -1, 0), "restricted_like_solution")
        self.assertNotEqual(m.classify_uhf(-1, -1, .9), "observed_lower_spin_broken_branch")
        self.assertNotEqual(m.classify_uhf(-1, -.8, 0), "observed_lower_spin_broken_branch")

    def test_planned_grid_is_unique_and_contains_prior_tests(self):
        grid = [r for _, r in m.point_plan()]
        self.assertEqual(len(grid), 25)
        self.assertEqual(len(set(grid)), 25)
        self.assertTrue(set([.62, .78, .94, 1.10, 1.26, 1.42, 1.58, 1.74, 1.9, 2.1, 2.3, 2.5, 2.7]).issubset(grid))

    def test_recovery_does_not_repeat_dispatched_or_successful_jobs(self):
        rows = [dict(basis="sto-3g", point_id=f"p{i}", method="FCI", status=status,
                     energy_driver_dispatched=called)
                for i, (status, called) in enumerate((("failed", False), ("failed", True), ("converged", True)))]
        self.assertEqual(m.recovery_targets(rows), {("sto-3g", "p0")})
        self.assertEqual(m.recovery_targets(rows[1:]), set())

    def test_fresh_study_aggregation_requires_no_recovery_directory(self):
        with tempfile.TemporaryDirectory(prefix="quantumequi-stages-") as directory:
            root = Path(directory)
            (root / "pilot").mkdir()
            for folder in (root, root / "pilot"):
                (folder / "energy_driver_jobs.json").write_text("[]", encoding="utf-8")
            self.assertEqual(m.study_stage_directories(root), [root / "pilot", root])
            (root / "pilot_fci_recovery").mkdir()
            (root / "pilot_fci_recovery/energy_driver_jobs.json").write_text("[]", encoding="utf-8")
            self.assertEqual(m.study_stage_directories(root),
                             [root / "pilot", root / "pilot_fci_recovery", root])

    def test_main_pilot_and_present_recovery_ledgers_are_mandatory(self):
        with tempfile.TemporaryDirectory(prefix="quantumequi-required-stages-") as directory:
            root = Path(directory)
            with self.assertRaises(FileNotFoundError):
                m.study_stage_directories(root)
            (root / "pilot").mkdir()
            (root / "pilot/energy_driver_jobs.json").write_text("[]", encoding="utf-8")
            with self.assertRaises(FileNotFoundError):
                m.study_stage_directories(root)
            (root / "energy_driver_jobs.json").write_text("[]", encoding="utf-8")
            (root / "pilot_fci_recovery").mkdir()
            with self.assertRaises(FileNotFoundError):
                m.study_stage_directories(root)


@unittest.skipUnless((OUT / "summary.json").exists(), "Run bounded study first.")
class SavedStudyChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
        cls.jobs = json.loads((OUT / "energy_driver_jobs.json").read_text(encoding="utf-8"))
        cls.curves = [j for j in cls.jobs if j["role"] == "curve"]
        cls.refs = json.loads((OUT / "equilibrium_reference.json").read_text(encoding="utf-8"))["references"]

    def test_exact_curve_and_roles(self):
        self.assertEqual(len(self.curves), 150)
        self.assertEqual(len({(j["basis"], j["method"], j["R_A"]) for j in self.curves}), 150)
        self.assertTrue(all(j["status"] == "converged" for j in self.jobs))
        self.assertEqual(sum(j["role"] == "isolated_H" for j in self.jobs), 2)
        self.assertEqual(sum(j["role"] == "curvature" for j in self.jobs), 10)

    def test_all_hashes_and_no_local_paths(self):
        for directory in (OUT, OUT / "pilot", OUT / "pilot_fci_recovery"):
            summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(m.sha(directory / summary["executed_source_snapshot"]), summary["execution"]["executed_source_sha256"])
            for relative, digest in summary["outputs_sha256"].items():
                self.assertEqual(m.sha(directory / relative), digest, str(directory / relative))
            for path in (directory / "logs").glob("*.out"):
                self.assertIsNone(re.search(r"[A-Za-z]:[\\/]", path.read_text(encoding="utf-8")))

    def test_driver_budget_and_pilot_failure_separation(self):
        all_jobs = list(self.jobs)
        for name in ("pilot", "pilot_fci_recovery"):
            all_jobs += json.loads((OUT / name / "energy_driver_jobs.json").read_text(encoding="utf-8"))
        self.assertLessEqual(sum(j["energy_driver_dispatched"] for j in all_jobs), 206)
        rejected = [j for j in all_jobs if not j["energy_driver_dispatched"]]
        self.assertEqual(len(rejected), 6)
        self.assertTrue(all(j["status"] == "failed" and "S_SQUARED" in j["error"] for j in rejected))
        self.assertEqual(self.summary["counts"]["energy_driver_calls"], len(self.jobs))

    def test_log_energy_and_scf_convergence_accounting(self):
        for job in self.jobs:
            text = (OUT / job["log_file"]).read_text(encoding="utf-8")
            self.assertEqual(text.count("Energy and wave function converged."), job["SCF_converged_messages"])
            pattern = r"Total CI energy\s*=\s*([-+0-9.]+)" if job["method"] == "FCI" else r"@(?:RHF|UHF) Final Energy:\s*([-+0-9.]+)"
            match = re.findall(pattern, text)
            self.assertTrue(match, job["job_id"])
            self.assertAlmostEqual(float(match[-1]), job["energy_Hartree"], places=10)

    def test_electronic_nuclear_total_identity(self):
        for job in self.jobs:
            self.assertAlmostEqual(job["energy_total_Hartree"], job["energy_electronic_Hartree"] + job["nuclear_repulsion_Hartree"], places=13)
            self.assertEqual(job["energy_total_Hartree"], job["energy_Hartree"])
            if job["R_A"] is not None:
                self.assertAlmostEqual(job["R_A"] * job["nuclear_repulsion_Hartree"], .52917721067, places=9)

    def test_variational_order_within_each_basis(self):
        for basis in m.BASES:
            for _, r in m.point_plan():
                e = {j["method"]: j["energy_Hartree"] for j in self.curves if j["basis"] == basis and j["R_A"] == r}
                self.assertLessEqual(e["FCI"], e["UHF"] + 1e-8)
                self.assertLessEqual(e["UHF"], e["RHF"] + 1e-8)

    def test_fci_determinant_count_from_log_not_broken_binding(self):
        for job in self.jobs:
            if job["method"] == "FCI":
                self.assertEqual(job["CI_determinants"], job["n_basis"] ** 2)
                self.assertIn(job["CI_determinants"], (4, 100))
                self.assertIsNone(job["S2"])
                self.assertIsNone(job["S2_source"])

    def test_reference_spin_recomputed_from_saved_orbitals(self):
        orbitals = json.loads((OUT / "reference_orbitals.json").read_text(encoding="utf-8"))
        by_id = {j["job_id"]: j for j in self.jobs}
        self.assertEqual(len(orbitals), len(self.jobs))
        for entry in orbitals:
            s2 = m.spin_squared_from_occupied(entry["occupied_alpha"], entry["occupied_beta"], entry["overlap"])
            self.assertAlmostEqual(s2, by_id[entry["job_id"]]["S2_reference_determinant"], places=12)

    def test_large_distance_uhf_branch_and_local_stability(self):
        for basis in m.BASES:
            group = {j["method"]: j for j in self.curves if j["basis"] == basis and j["R_A"] == 4}
            self.assertGreater(group["UHF"]["S2"], .99)
            self.assertLess(group["UHF"]["energy_Hartree"], group["RHF"]["energy_Hartree"] - .05)
            self.assertTrue(group["UHF"]["stability_eigenvalues"])
            self.assertGreaterEqual(min(group["UHF"]["stability_eigenvalues"]["SCF STABILITY EIGENVALUES"]), -1e-7)

    def test_atomic_doublets_and_dissociation_reference(self):
        for ref in self.refs:
            atom = next(j for j in self.jobs if j["role"] == "isolated_H" and j["basis"] == ref["basis"])
            self.assertEqual(atom["multiplicity_input"], 2)
            self.assertAlmostEqual(atom["S2"], .75)
            self.assertEqual(ref["E_infinity_Hartree"], 2 * atom["energy_Hartree"])
            self.assertAlmostEqual(ref["D_e_Hartree"], ref["E_infinity_Hartree"] - ref["E_min_Hartree"], places=14)

    def test_minimum_curvature_independent_stencil(self):
        by_id = {j["job_id"]: j for j in self.jobs}
        for ref in self.refs:
            self.assertTrue(ref["optimization_success"], ref["optimization_message"])
            self.assertGreater(ref["R_e_A"], .6)
            self.assertLess(ref["R_e_A"], .9)
            local = [by_id[j] for j in ref["curvature_job_ids"]]
            e = np.array([j["energy_Hartree"] for j in local])
            h = ref["curvature_step_A"]
            expected = np.dot(e, [-1, 16, -30, 16, -1]) / (12*h*h)
            self.assertAlmostEqual(ref["curvature_Hartree_A2"], expected, places=9)
            self.assertGreater(expected, 0)
            self.assertLess(abs(ref["gradient_5point_Hartree_A"]), 1e-5)

    def test_curve_csv_is_lossless_public_interface(self):
        with (OUT / "curve.csv").open(encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))
        by_id = {j["job_id"]: j for j in self.curves}
        self.assertEqual(len(rows), 150)
        for row in rows:
            source = by_id[row["job_id"]]
            self.assertEqual(float(row["energy_Hartree"]), source["energy_Hartree"])
            self.assertEqual(float(row["R_A"]), source["R_A"])
            self.assertEqual(row["status"], "converged")


if __name__ == "__main__":
    unittest.main()
