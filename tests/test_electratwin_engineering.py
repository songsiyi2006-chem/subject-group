from pathlib import Path
import sys,unittest
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'electratwin/scripts'))
import engineering_audit as e

class ElectraTwinEngineeringTests(unittest.TestCase):
    def test_geometry_flow_and_poiseuille_integral(self):
        r=e.hydraulic();self.assertAlmostEqual(r['volume_mL'],.216);self.assertAlmostEqual(r['residence_time_s'],28.8)
        pressure_gradient=r['parallel_plate_pressure_drop_Pa']/r['length_m'];mu=r['viscosity_Pa_s'];H=r['height_m'];W=r['width_m']
        q=W*quad(lambda y:pressure_gradient/(2*mu)*y*(H-y),0,H)[0]
        self.assertAlmostEqual(q/r['flow_m3_s'],1,places=10)

    def test_gap_scaling_and_reynolds_units(self):
        a=e.hydraulic();b=e.hydraulic(height_m=.0006)
        self.assertAlmostEqual(a['parallel_plate_pressure_drop_Pa']/b['parallel_plate_pressure_drop_Pa'],8)
        self.assertAlmostEqual(b['residence_time_s']/a['residence_time_s'],2)
        self.assertGreater(a['Re_hydraulic'],2);self.assertLess(a['Re_hydraulic'],3)

    def test_heat_capacity_balance_and_invalid_inputs(self):
        delta=e.electrical_heating(.1,3,1e-8,density_kg_m3=1000,heat_capacity_J_kg_K=3000)
        self.assertAlmostEqual(delta,10);self.assertEqual(e.electrical_heating(0,3,1e-8),0)
        for kwargs in [{'flow_uL_min':0},{'height_m':-1},{'viscosity_Pa_s':float('nan')}]:
            with self.assertRaises(ValueError):e.hydraulic(**kwargs)
        with self.assertRaises(ValueError):e.electrical_heating(.1,3,1e-8,heat_fraction=2)

if __name__=='__main__':unittest.main()
