import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('error_budget',Path(__file__).resolve().parents[1]/'quantumequi/scripts/error_budget_extension.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class ErrorBudgetTest(unittest.TestCase):
    def test_signed_cancellation_is_not_zero_components(self):
        r=module.energy_components(-1.1,-1.0,-1.0,-1.1)
        self.assertEqual(r['total_error_vs_FCI_Hartree'],0)
        self.assertLess(r['learning_error_Hartree'],0)
        self.assertGreater(r['RHF_minus_FCI_Hartree'],0)
        self.assertAlmostEqual(r['reconstruction_residual_Hartree'],0)

    def test_separate_reference_replay_drift(self):
        r=module.energy_components(-1.05,-1.0,-1.001,-1.1)
        self.assertAlmostEqual(r['reference_replay_drift_Hartree'],.001)
        self.assertAlmostEqual(r['reconstruction_residual_Hartree'],0)

    def test_common_energy_zero_does_not_change_components(self):
        a=module.energy_components(.1,.2,.21,.05)
        b=module.energy_components(10.1,10.2,10.21,10.05)
        for key in a:self.assertAlmostEqual(a[key],b[key],places=13)

    def test_invalid_reference_rejected(self):
        with self.assertRaises(ValueError):module.energy_components(1,float('nan'),2,3)


if __name__=='__main__':unittest.main()
