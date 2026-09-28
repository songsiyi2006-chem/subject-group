"""Symmetry, differentiation and structural checks on the archived source model."""
import json
from pathlib import Path
import unittest
import numpy as np
import torch
from quantumequi.scripts import potential_audit as audit


class PotentialAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)
        captured = json.loads((audit.BASE/'results/original/captured_results.json').read_text())
        cls.model = audit.load_model(torch.float64)
        cls.z = np.array(captured['z_atomic_numbers'])
        cls.x = np.array(captured['R_reactant'])

    def test_rigid_covariance_including_reflection(self):
        e, f = audit.evaluate(self.model, self.z, self.x)
        q = audit.orthogonal(199, reflection=True)
        e2, f2 = audit.evaluate(self.model, self.z, self.x@q.T+[1.2,-.4,2.5])
        self.assertLess(abs(e2-e), 1e-10)
        np.testing.assert_allclose(f2, f@q.T, atol=1e-11, rtol=1e-10)

    def test_zero_total_force_and_torque(self):
        _, f = audit.evaluate(self.model, self.z, self.x)
        np.testing.assert_allclose(f.sum(0), 0, atol=1e-11)
        np.testing.assert_allclose(np.cross(self.x, f).sum(0), 0, atol=1e-11)

    def test_force_matches_central_difference(self):
        _, force = audit.evaluate(self.model, self.z, self.x)
        fd = audit.finite_difference(self.model, self.z, self.x, 1e-4)
        np.testing.assert_allclose(fd, force, atol=1e-8, rtol=1e-6)

    def test_atomic_relabeling_permutes_force(self):
        order = np.array([7,4,2,5,0,6,1,3])
        e, f = audit.evaluate(self.model, self.z, self.x)
        e2, f2 = audit.evaluate(self.model, self.z[order], self.x[order])
        self.assertAlmostEqual(e2, e, places=11)
        np.testing.assert_allclose(f2, f[order], atol=1e-11, rtol=1e-10)


if __name__ == '__main__':
    unittest.main()
