"""Analytic variance decompositions, sampling contracts and bounded PDE checks."""
import unittest
from unittest.mock import patch

import numpy as np

from electratwin.scripts.uncertainty_analysis import (
    PARAMETERS, OUTPUTS, transform_unit, make_design, jansen_indices,
    paired_bootstrap, metrics_from_summary, evaluate_point, local_design,
)


class TestElectraTwinUncertainty(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a, cls.b, cls.ab = make_design(power=12, seed=71)

    def test_bounds_and_log_uniform_midpoint(self):
        x = transform_unit(np.array([[0.]*5, [.5]*5, [1.]*5]))
        np.testing.assert_allclose(x[0], [p["low"] for p in PARAMETERS], rtol=1e-14)
        np.testing.assert_allclose(x[2], [p["high"] for p in PARAMETERS], rtol=1e-14)
        self.assertAlmostEqual(x[1,3], .08)
        np.testing.assert_allclose(x[1,[0,1,2,4]], [450., .48, 1.1e-9, .0003])

    def test_invalid_unit_coordinates_rejected(self):
        for x in [[.5]*4, [.5]*4+[float("nan")], [.5]*4+[-.001], [.5]*4+[1.001]]:
            with self.subTest(x=x), self.assertRaises(ValueError):
                transform_unit(x)

    def test_scrambled_design_reproducible_and_prefix_preserved(self):
        a,b,ab=make_design(4,381)
        a2,b2,ab2=make_design(5,381)
        np.testing.assert_array_equal(a,a2[:16]); np.testing.assert_array_equal(b,b2[:16])
        np.testing.assert_array_equal(ab,ab2[:,:16])
        self.assertTrue(np.all((a>=0)&(a<1)))
        self.assertFalse(np.array_equal(a,make_design(4,382)[0]))

    def test_pick_freeze_replaces_exactly_one_column(self):
        for i in range(5):
            np.testing.assert_array_equal(self.ab[i,:,i],self.b[:,i])
            keep=[j for j in range(5) if j!=i]
            np.testing.assert_array_equal(self.ab[i][:,keep],self.a[:,keep])

    def test_additive_analytic_sobol_indices(self):
        coefficients=np.array([1.,2.,3.,0.,0.])
        result=jansen_indices(self.a@coefficients,self.b@coefficients,self.ab@coefficients)
        expected=coefficients**2/14
        np.testing.assert_allclose(result["first_order"][:,0],expected,atol=.005)
        np.testing.assert_allclose(result["total_order"][:,0],expected,atol=.005)

    def test_product_analytic_interaction_indices(self):
        f=lambda x:x[...,0]*x[...,1]
        result=jansen_indices(f(self.a),f(self.b),f(self.ab))
        np.testing.assert_allclose(result["first_order"][:,0],[3/7,3/7,0,0,0],atol=.006)
        np.testing.assert_allclose(result["total_order"][:,0],[4/7,4/7,0,0,0],atol=.006)

    def test_constant_output_has_no_fabricated_indices(self):
        with self.assertRaises(ValueError):
            jansen_indices(np.ones(8),np.ones(8),np.ones((5,8)))

    def test_finite_sample_negative_and_above_one_indices_retained(self):
        a=np.array([0.,1.,0.,1.]); b=np.array([1.,0.,1.,0.]); ab=np.full((1,4),10.)
        result=jansen_indices(a,b,ab)
        self.assertLess(result["first_order"][0,0],0)
        self.assertGreater(result["total_order"][0,0],1)

    def test_bootstrap_resamples_the_same_rows_in_all_blocks(self):
        a,b,ab=make_design(4,812)
        f=lambda x:x[...,0]+2*x[...,2]
        ya,yb,yab=f(a),f(b),f(ab)
        result=paired_bootstrap(ya,yb,yab,replicates=3,seed=91)
        rows=np.random.default_rng(91).integers(0,len(ya),size=len(ya))
        expected=jansen_indices(ya[rows],yb[rows],yab[:,rows])
        np.testing.assert_array_equal(result["first_order"][0],expected["first_order"])
        np.testing.assert_array_equal(result["total_order"][0],expected["total_order"])

    def test_productivity_and_energy_units(self):
        r=1e-7; current=2*96485.33212*r
        summary=dict(converged=True,nonnegative=True,linear_residual_relative=0.,
            material_balance_relative_error=0.,charge_balance_relative_error=0.,
            reacted_mol_s=r,current_A=current,reactor_volume_m3=1e-6,conversion_pct=50.)
        metrics=metrics_from_summary(summary,.48)
        kg_day=r*331.43/1000*60*60*24
        voltage=1.85+.48+18.5*current
        kw=voltage*current/1000; kg_h=r*331.43/1000*60*60
        self.assertAlmostEqual(metrics["STY_assumed_product_kg_m3_day"],kg_day/1e-6,places=9)
        self.assertAlmostEqual(metrics["SEC_assumed_product_kwh_kg"],kw/kg_h,places=12)

    def test_bounded_PDE_outputs_are_positive_and_conservative(self):
        for unit in [np.full(5,.5),np.array([0,1,0,1,0])]:
            row=evaluate_point(unit,nx=12,ny=8)
            self.assertTrue(row["converged"])
            self.assertTrue(all(row[name]>0 for name in OUTPUTS))
            self.assertLess(row["material_balance_relative_error"],1e-10)
            self.assertLess(row["charge_balance_relative_error"],1e-10)

    def test_local_design_has_twenty_one_symmetric_quantile_points(self):
        points=local_design()
        self.assertEqual(len(points),21)
        self.assertEqual(points[0][0],"center")
        for index in range(1,len(points),2):
            left,right=points[index][-1],points[index+1][-1]
            np.testing.assert_allclose((left+right)/2,.5)
            self.assertTrue(np.all((left>=0)&(left<=1)))

    def test_failed_PDE_is_rejected_before_statistics(self):
        with patch("electratwin.scripts.uncertainty_analysis.solve_transport",return_value={"summary":{"converged":False,"nonnegative":True}}):
            with self.assertRaises(ValueError):
                evaluate_point(np.full(5,.5))


if __name__=="__main__":
    unittest.main()
