"""Independent analytical invariants and provenance gates; no assay certification."""
from pathlib import Path
import csv,hashlib,json,sys,tempfile,unittest
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];TK=ROOT/'toolkit';sys.path.insert(0,str(TK/'scripts'))
import audit_toolkit as a

def rows(name):
    with (TK/name).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

class ToolkitTests(unittest.TestCase):
    def test_hplc_units_from_known_amount(self):
        # 0.1 mmol in 10 mL = 10 mM; after 5-fold dilution, 2 mM.
        r=a.hplc(2007,1000,7,5,10,.2,(0,3))
        self.assertAlmostEqual(r['product_mmol'],.1);self.assertAlmostEqual(r['yield_pct'],50)
        self.assertTrue(r['numeric_range_pass']);self.assertFalse(r['measurement_validated'])

    def test_source_abnormal_yield_is_never_clipped(self):
        r=a.hplc(245600,14250,120,25,6,.2)
        self.assertAlmostEqual(r['yield_pct'],1292);self.assertEqual(len(r['flags']),2)
        self.assertFalse(r['numeric_range_pass']);self.assertGreater(r['product_mmol'],.2)

    def test_negative_and_out_of_range_signal_flags(self):
        self.assertIn('below_intercept_not_a_negative_product_amount',a.hplc(0,1000,1,1,1,1)['flags'])
        self.assertIn('outside_calibration_range',a.hplc(4000,1000,0,1,1,1,(0,2))['flags'])

    def test_invalid_analytical_denominators(self):
        for value in [0,-1,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):a.hplc(100,value,0,1,1,1)
            with self.assertRaises(ValueError):a.qnmr(1,value,1,1,1,100,1,1)
        for purity in [0,1.1,float('nan')]:
            with self.assertRaises(ValueError):a.qnmr(1,1,1,1,1,100,purity,1)

    def test_nmr_proton_mass_purity_and_aliquot_units(self):
        # Integral-per-proton ratio = 2; standard = 0.1 mmol at 100% purity.
        base=a.qnmr(4,3,2,3,10,100,1,1)
        self.assertAlmostEqual(base['product_mmol'],.2)
        self.assertAlmostEqual(a.qnmr(4,3,2,3,10,100,.9,1,.5)['product_mmol'],.36)
        self.assertEqual(a.qnmr(0,3,2,3,10,100,1,1)['yield_pct'],0)
        self.assertFalse(base['measurement_validated'])

    def test_nmr_rounding_separate_from_hplc_discrepancy(self):
        exact=a.qnmr(5.12,3,1,3,6,168.19,1,.2)['yield_pct']
        rounded=5.12*.0357/.2*100
        self.assertAlmostEqual(exact,91.32528687793567)
        self.assertLess(abs(exact-rounded),.1);self.assertGreater(1292-rounded,1200)

    def test_roles_units_ids_and_partial_range_rejected(self):
        source=rows('data/source_integrated_area.csv')[0]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'data.csv'
            for change in [{'evidence_role':'experimental'},{'volume_unit':'L'},{'area_mode':'ratio'},{'sample_id':''},{'calibration_min_mM':'0'},{'calibration_max_mM':'2'}]:
                a.cwrite(p,[dict(source,**change)])
                with self.assertRaises(ValueError):a.parse_integrated_csv(p,'source_example')
            a.cwrite(p,[source,source])
            with self.assertRaises(ValueError):a.parse_integrated_csv(p,'source_example')

    def experimental_hplc(self,folder):
        source=rows('data/source_integrated_area.csv')[0]
        row=dict(source,evidence_role='experimental',product_identity='test-only identity',analyst='test fixture',recorded_at='2026-09-28T00:00:00Z',calibration_min_mM='0',calibration_max_mM='20')
        for kind in ['raw','calibration']:
            p=folder/(kind+'.txt');p.write_text('Unit-test placeholder, not real analytical evidence.')
            row[kind+'_file']=p.name;row[kind+'_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
        return row

    def test_experimental_hash_gate_does_not_certify_measurement(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);row=self.experimental_hplc(folder);p=folder/'records.csv';a.cwrite(p,[row])
            result=a.parse_integrated_csv(p,'experimental')[0]
            self.assertFalse(result['measurement_validated']);self.assertFalse(result['numeric_range_pass'])
            (folder/'raw.txt').write_text('Changed bytes.')
            with self.assertRaises(ValueError):a.parse_integrated_csv(p,'experimental')

    def test_experimental_path_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent=Path(tmp);folder=parent/'data';folder.mkdir();row=self.experimental_hplc(folder)
            outside=parent/'outside.txt';outside.write_text('Outside evidence.')
            row.update(raw_file='../outside.txt',raw_sha256=hashlib.sha256(outside.read_bytes()).hexdigest())
            p=folder/'records.csv';a.cwrite(p,[row])
            with self.assertRaises(ValueError):a.parse_integrated_csv(p,'experimental')

    def test_nmr_metadata_and_certificate_required(self):
        source=rows('data/source_nmr_integrals.csv')[0]
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);p=folder/'nmr.csv';r=dict(source,evidence_role='experimental');a.cwrite(p,[r])
            with self.assertRaises(ValueError):a.parse_nmr_csv(p,'experimental')
            r.update(product_identity='test-only',product_peak_assignment='test-only',relaxation_record='test-only',analyst='fixture',recorded_at='2026-09-28T00:00:00Z')
            for kind in ['raw','standard_certificate']:
                q=folder/(kind+'.txt');q.write_text('Test placeholder, not scientific evidence.');r[kind+'_file']=q.name;r[kind+'_sha256']=hashlib.sha256(q.read_bytes()).hexdigest()
            a.cwrite(p,[r]);self.assertFalse(a.parse_nmr_csv(p,'experimental')[0]['measurement_validated'])
            r['standard_certificate_sha256']='0'*64;a.cwrite(p,[r])
            with self.assertRaises(ValueError):a.parse_nmr_csv(p,'experimental')

    def test_nmr_source_units_and_identity(self):
        source=rows('data/source_nmr_integrals.csv')[0]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'data.csv'
            for change in [{'mass_unit':'g'},{'amount_unit':'mol'},{'evidence_role':'synthetic'},{'product_integral':'-1'},{'sample_fraction':'0'}]:
                a.cwrite(p,[dict(source,**change)])
                with self.assertRaises(ValueError):a.parse_nmr_csv(p,'source_example')
            a.cwrite(p,[source,source])
            with self.assertRaises(ValueError):a.parse_nmr_csv(p,'source_example')

    def test_gaussian_area_has_independent_quadrature_normalization(self):
        integral=quad(lambda t:a.gaussian(t,6.42,.07),5,8,epsabs=1e-10)[0]
        self.assertAlmostEqual(integral,1,places=9)

    def test_known_shape_noiseless_area_recovery(self):
        t=np.linspace(5.5,7.5,2001);signal=50+8*(t-6.5)+1234*a.gaussian(t,6.42,.07)+678*a.gaussian(t,6.57,.09)
        r=a.deconvolve_known_shapes(t,signal,6.42,6.57,.07,.09)
        self.assertAlmostEqual(r['area1'],1234,places=8);self.assertAlmostEqual(r['area2'],678,places=8)

    def test_peak_width_negative_result_reproduces(self):
        t=np.linspace(5.8,7.2,1401);signal=50+8*(t-6.5)+17500*a.gaussian(t,6.42,.08)+7000*a.gaussian(t,6.45,.09)+np.random.default_rng(0).normal(0,300,len(t))
        r=a.deconvolve_known_shapes(t,signal,6.42,6.45,.07,.09)
        self.assertGreater(abs(r['area1']/17500-1),.2)
        saved=next(r for r in rows('results/audit/deconvolution_recovery.csv') if r['width_mismatch']=='True' and r['seed']=='0' and r['separation_min']=='0.03')
        self.assertAlmostEqual(r['area1'],float(saved['fitted_area1']),places=7)

    def test_pareto_ties_dominance_and_nonfinite(self):
        values=[[1,1],[1,1],[2,2],[0,3],[3,0]]
        np.testing.assert_equal(a.pareto_mask(values),[True,True,False,True,True])
        with self.assertRaises(ValueError):a.pareto_mask([[1,np.nan]])

    def test_pareto_saved_front_matches_independent_loop(self):
        r=rows('results/audit/pareto_conditions.csv');v=[[-float(x['yield_pct']),-float(x['FE_pct']),float(x['SEC_kWh_kg'])] for x in r]
        expected=[]
        for p in v:
            dominated=any(all(qi<=pi for qi,pi in zip(q,p)) and any(qi<pi for qi,pi in zip(q,p)) for q in v)
            expected.append(not dominated)
        self.assertEqual(sum(expected),23);self.assertEqual(expected,[x['nondominated']=='True' for x in r])
        self.assertEqual(sum(x['source_threshold']=='True' for x in r),0)

    def test_energy_from_independent_charge_mass_bookkeeping(self):
        amount_mol=.002;charge=2*96485.33*amount_mol/.8;energy_kWh=3.2*charge/3.6e6;mass_kg=amount_mol*.223
        d=json.loads((TK/'results/audit/figure_audit.json').read_text())
        self.assertAlmostEqual(energy_kWh/mass_kg,d['Pareto']['energy_numerator_kWh_kg']/.8)
        self.assertFalse(d['energy']['DFT_executed']);self.assertFalse(d['scope']['chemical_prediction_validated'])

    def test_synthetic_calibration_fit_has_normal_equation_residuals(self):
        r=rows('results/audit/synthetic_calibration.csv');x=np.array([float(v['concentration_mM']) for v in r]);y=np.array([float(v['area']) for v in r]);d=json.loads((TK/'results/audit/synthetic_recovery.json').read_text())
        slope=np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2)
        self.assertAlmostEqual(slope,d['calibration_slope'],places=8)
        self.assertAlmostEqual(y.mean()-slope*x.mean(),d['calibration_intercept'],places=8)
        self.assertTrue(all(v['evidence_role']=='synthetic' for v in r))

    def test_budget_and_unsupported_grant_status(self):
        r=rows('results/audit/grant_budget.csv');total=sum(float(v['quantity'])*float(v['assumed_unit_CNY']) for v in r)
        self.assertEqual(total,10000)
        d=json.loads((TK/'results/audit/grant_planning.json').read_text())
        self.assertFalse(d['approval_or_grant_submission']);self.assertIsNone(d['confirmed_funding_CNY'])

if __name__=='__main__':unittest.main()
