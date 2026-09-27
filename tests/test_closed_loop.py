"""Independent numerical and evidence-boundary tests for the closed-loop extension."""
from pathlib import Path
import copy,csv,hashlib,json,sys,tempfile,unittest
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];CL=ROOT/'closed_loop';sys.path.insert(0,str(CL/'scripts'))
import audit_closed_loop as a
import generate_reviewed_inputs as g

def data(p):return json.loads((CL/p).read_text())
def rows(p):
    with (CL/p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

class ClosedLoopTests(unittest.TestCase):
    def test_charge_units_and_independent_integral(self):
        s=a.stoichiometry();self.assertAlmostEqual(s['charge_C'],42.4535452,places=9)
        self.assertAlmostEqual(quad(lambda t:.01875,0,s['duration_min']*60)[0],s['charge_C'],places=10)
        self.assertAlmostEqual(s['substrate_mg'],44.655,places=6);self.assertAlmostEqual(s['electrolyte_mmol'],.6)

    def test_parameter_changes_propagate(self):
        base=a.stoichiometry();doubleV=a.stoichiometry(volume_mL=12);doubleA=a.stoichiometry(area_cm2=3);doubleN=a.stoichiometry(scale_mmol=.4)
        self.assertAlmostEqual(doubleV['MeCN_mL'],9.6);self.assertAlmostEqual(doubleV['HFIP_mL'],2.4)
        self.assertAlmostEqual(doubleV['substrate_M'],base['substrate_M']/2);self.assertAlmostEqual(doubleV['electrolyte_mg'],2*base['electrolyte_mg'])
        self.assertAlmostEqual(doubleA['duration_min'],base['duration_min']/2);self.assertAlmostEqual(doubleN['duration_min'],2*base['duration_min'])

    def test_invalid_stoichiometry_rejected(self):
        for kwargs in [{'scale_mmol':0},{'area_cm2':-1},{'j_mA_cm2':float('nan')},{'meCN_fraction':1}]:
            with self.assertRaises(ValueError):a.stoichiometry(**kwargs)

    def test_two_electron_material_balance(self):
        s=a.stoichiometry();self.assertAlmostEqual(s['full_conversion_FE_ceiling_pct'],100*2/2.2)
        for r in rows('results/audit/charge_balance_scenarios.csv'):
            self.assertAlmostEqual(float(r['FE_pct'])/100*s['charge_C']/(2*a.F)*1000,float(r['product_mmol']))

    def test_partner_identity_not_interchangeable(self):
        r=rows('results/audit/partner_scenarios.csv');one=[v for v in r if float(v['assumed_equivalents'])==1]
        self.assertEqual({v['formula'] for v in one},{'C7H8S','C6H6S'})
        self.assertNotEqual(one[0]['mass_mg'],one[1]['mass_mg']);self.assertFalse(data('results/audit/stoichiometry.json')['ready_for_wet_lab'])

    def test_given_mae_is_reproduced_from_source_arrays(self):
        X,p,y=a.source_arrays();d=data('results/audit/active_learning.json');self.assertEqual(X.shape,(6,3));self.assertAlmostEqual(np.abs(p-y).mean(),17.3/6);self.assertAlmostEqual(d['given_prediction_MAE_pp'],17.3/6)

    def test_loo_mean_baseline_independent(self):
        X,p,y=a.source_arrays();expected=np.array([(sum(y)-v)/5 for v in y]);self.assertAlmostEqual(np.abs(expected-y).mean(),13.98)
        r=[v for v in rows('results/audit/loo_predictions.csv') if v['method']=='mean_baseline'];np.testing.assert_allclose([float(v['predicted_yield_pct']) for v in r],expected)

    def test_one_heldout_gp_fit_reproduces_saved_prediction(self):
        X,p,y=a.source_arrays();gp,sc,_=a.fit_gp(X[1:],y[1:],'scaled_mle');mu,_,_=a.predictions(gp,sc,X[:1])
        r=next(v for v in rows('results/audit/loo_predictions.csv') if v['method']=='scaled_mle' and v['held_out_run']=='mock-01');self.assertAlmostEqual(mu[0],float(r['predicted_yield_pct']),places=6)
        np.testing.assert_allclose(sc.mean_,X[1:].mean(axis=0))

    def test_ei_matches_independent_quadrature(self):
        for mu,sigma,best in [(78.,12.,88.),(91.,2.,88.),(88.,5.,88.)]:
            direct=quad(lambda x:max(x-best,0)*np.exp(-.5*((x-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi)),best,mu+12*sigma,epsabs=1e-10)[0]
            self.assertAlmostEqual(float(a.expected_improvement(mu,sigma,best)),direct,places=8)
        self.assertEqual(float(a.expected_improvement(90,0,88)),2);self.assertEqual(float(a.expected_improvement(80,0,88)),0)

    def test_noise_variance_removed_in_physical_units(self):
        X,p,y=a.source_arrays();gp,sc,_=a.fit_gp(X,y,'scaled_fixed',2);mu,latent,obs=a.predictions(gp,sc,a.candidate_pool());np.testing.assert_allclose(obs**2-latent**2,4,atol=1e-10)

    def test_candidate_repeats_excluded(self):
        X,_,_=a.source_arrays();pool=a.candidate_pool();mask=a.unseen_mask(pool,X);self.assertEqual(len(pool),90);self.assertEqual(sum(mask),89)
        d=data('results/audit/active_learning.json')
        for methods in d['selections'].values():
            for ranked in methods.values():
                for r in ranked:self.assertTrue(mask[r['candidate']])

    def test_negative_model_comparison_retained(self):
        d=data('results/audit/active_learning.json');self.assertGreater(d['LOO']['scaled_mle']['MAE_pp'],d['LOO']['mean_baseline']['MAE_pp'])
        self.assertEqual(d['evidence_role'],'mock');mc=d['EI_monte_carlo'];self.assertLess(abs(mc['estimated_pp']-mc['analytic_pp']),4*mc['standard_error_pp'])

    def test_evidence_roles_and_units_reject_invalid_rows(self):
        original=rows('data/mock_feedback.csv')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'records.csv'
            for key,value in [('yield_pct','101'),('temperature_K','30'),('yield_pct','nan'),('run_id','mock-02'),('reaction_id','other')]:
                bad=copy.deepcopy(original);bad[0][key]=value;a.save_csv(path,bad)
                with self.assertRaises(ValueError):a.load_feedback(path)
            a.save_csv(path,original)
            with self.assertRaises(ValueError):a.load_feedback(path,'experimental')

    def test_experimental_linkage_and_hash_guard(self):
        # Test fixtures are not measurements and are confined to a temporary directory.
        original=rows('data/mock_feedback.csv')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);raw=root/'fixture.dat';raw.write_bytes(b'test fixture, not laboratory evidence');digest=hashlib.sha256(raw.read_bytes()).hexdigest()
            for r in original:
                r.update(evidence_role='experimental',reaction_id='test_only',substrate_smiles='CCO',partner_smiles='O',product_smiles='CC=O',assay='HPLC',operator='test-fixture',recorded_at='test-time')
                for kind in ['assay','calibration','electrolysis']:r[kind+'_file']='fixture.dat';r[kind+'_sha256']=digest
            path=root/'records.csv';a.save_csv(path,original);self.assertEqual(len(a.load_feedback(path,'experimental')),6)
            raw.write_bytes(b'changed')
            with self.assertRaises(ValueError):a.load_feedback(path,'experimental')
            raw.write_bytes(b'test fixture, not laboratory evidence');original[0]['assay_file']='../fixture.dat';a.save_csv(path,original)
            with self.assertRaises(ValueError):a.load_feedback(path,'experimental')

    def test_quantum_identity_spin_and_geometry(self):
        manifest=data('inputs/reviewed/input_manifest.json');self.assertEqual(len(manifest['files']),12);self.assertFalse(manifest['dft_executed'])
        self.assertEqual(manifest['molecules']['target_indole']['neutral_electrons'],118);self.assertEqual(manifest['molecules']['indoline_control']['neutral_electrons'],64)
        for smiles,folder,identity in g.MOLECULES.values():
            coords,info=g.read_verified_geometry(smiles,folder,identity);self.assertEqual(len(coords.splitlines()),info['atoms']);self.assertGreater(info['minimum_pair_distance_A'],.6)
        for r in manifest['files']:
            self.assertEqual((r['electrons']-(r['multiplicity']-1))%2,0)
            p=CL/'inputs/reviewed'/r['path'];self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r['sha256'])

    def test_dft_solvation_and_dependent_geometry(self):
        for p in (CL/'inputs/reviewed').glob('*.inp'):
            text=p.read_text();self.assertIn('smd true',text);self.assertIn('SMDsolvent "acetonitrile"',text)
            if p.stem.endswith('_sp'):
                self.assertIn('M062X def2-TZVP',text);self.assertIn('* xyzfile',text);self.assertNotIn('D3BJ',text)
                required=text.split('* xyzfile')[1].splitlines()[0].split()[-1];self.assertFalse((p.parent/required).exists())
            else:self.assertIn('B3LYP/G D3BJ',text)
        for p in (CL/'inputs/reviewed').glob('*.gjf'):
            text=p.read_text();self.assertEqual(text.count('--Link1--'),1);self.assertEqual(text.count('solvent=acetonitrile'),2);self.assertIn('%mem=3GB',text)

    def test_original_source_and_execution_hash(self):
        source=data('source/source_record.json');run=data('results/original/execution.json')
        self.assertEqual(hashlib.sha256((CL/'source/run_closed_loop_platform.py').read_bytes()).hexdigest(),source['script_sha256']);self.assertEqual(run['script_sha256'],source['script_sha256']);self.assertEqual(run['returncode'],0)
        self.assertFalse(run['dft_executed']);self.assertFalse(run['wet_lab_executed'])

if __name__=='__main__':unittest.main()
