"""Scientific invariants and source-to-result checks; no xTB rerun required."""
import hashlib,json,sys,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score,mean_squared_error
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from extended_benchmarks import sequential_solution,flow_metrics,redox_data,classify_label

class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ext=ROOT/'results/extended';cls.mol=ROOT/'results/molecular'
        cls.summary=json.loads((cls.ext/'summary.json').read_text())

    def test_original_bytes_preserved(self):
        for name,expected in [('simulate_all_topics.py','41caab8993548814500ff447aa881bf04e62581a418f160727d7e78a3e1ef389'),('benchmark_results.json','1775377f34aecba1dbb5fe17fe0b9a6625b6de3a81a31f66247e694317f2ea41')]:
            self.assertEqual(hashlib.sha256((ROOT/'data/original'/name).read_bytes()).hexdigest(),expected)

    def test_bo_matched_design_and_budget(self):
        df=pd.read_csv(self.ext/'bo_observations.csv')
        self.assertEqual(len(df),30*3*20)
        for seed,g in df.groupby('seed'):
            initial=[v.head(8)[['j','concentration','epsilon','observed']].to_numpy() for _,v in g.groupby('method')]
            np.testing.assert_array_equal(initial[0],initial[1]);np.testing.assert_array_equal(initial[1],initial[2])
            self.assertTrue((g.groupby('method').size()==20).all())

    def test_bo_regret_computed_from_recommended_point(self):
        obs=pd.read_csv(self.ext/'bo_observations.csv');runs=pd.read_csv(self.ext/'bo_runs.csv')
        for row in runs.itertuples():
            g=obs[(obs.seed==row.seed)&(obs.method==row.method)]
            rec=g.loc[g.observed.idxmax()]
            self.assertAlmostEqual(92-rec.latent,row.recommendation_regret,places=10)
            self.assertGreaterEqual(row.recommendation_regret+1e-12,row.oracle_simple_regret)

    def test_bo_pairwise_interval_data(self):
        from extended_benchmarks import paired_bootstrap
        p=pd.read_csv(self.ext/'bo_runs.csv').pivot(index='seed',columns='method',values='recommendation_regret')
        for m in ['original_gp','scaled_gp']:
            r=paired_bootstrap(p[m]-p.random)
            self.assertEqual(r,self.summary['bo']['paired_recommendation_regret'][m+'_minus_random'])

    def test_regression_test_not_in_train(self):
        test=pd.read_csv(self.ext/'regression_test.csv').iloc[:,:4].to_numpy()
        testset={tuple(x) for x in np.round(test,12)}
        for seed in range(30):
            for n in [20,60,120,240]:
                train,_=redox_data(n,np.random.default_rng(40000+seed*1000+n))
                self.assertFalse(testset.intersection(tuple(x) for x in np.round(train,12)))

    def test_regression_independent_refit(self):
        test=pd.read_csv(self.ext/'regression_test.csv');x,y=redox_data(60,np.random.default_rng(40060))
        pred=LinearRegression().fit(x,y).predict(test.iloc[:,:4].to_numpy())
        rmse=np.sqrt(mean_squared_error(test.synthetic_target,pred))
        rows=pd.read_csv(self.ext/'regression_runs.csv')
        stored=rows[(rows.seed==0)&(rows.n_train==60)&(rows.model=='linear')].test_rmse.iloc[0]
        self.assertAlmostEqual(rmse,stored,places=10)

    def test_sac_probability_partition_and_symmetry(self):
        df=pd.read_csv(self.ext/'sac_rank_sensitivity.csv')
        for _,g in df.groupby('noise_sd'):
            self.assertAlmostEqual(g.probability_rank1.sum(),1,places=12)
        g=df[df.noise_sd==.4].set_index('site')
        for c in ['N4','N3C1','N2O2']:
            self.assertEqual(g.loc['Co-'+c].noiseless_mean,g.loc['Cu-'+c].noiseless_mean)
            self.assertLess(abs(g.loc['Co-'+c].probability_rank1-g.loc['Cu-'+c].probability_rank1),.01)

    def test_classification_refit_and_metrics(self):
        x=np.random.default_rng(70000).uniform(.1,2.5,(80,4));y=classify_label(x)
        test=np.random.default_rng(777).uniform(.1,2.5,(5000,4));truth=classify_label(test)
        p=RandomForestClassifier(n_estimators=30,random_state=42).fit(x,y).predict_proba(test)[:,1]
        saved=pd.read_csv(self.ext/'classification_test_predictions_seed0.csv')
        np.testing.assert_array_equal(truth,saved.y);np.testing.assert_allclose(p,saved.probability,atol=1e-14)
        rows=pd.read_csv(self.ext/'classification_runs.csv');v=rows[(rows.seed==0)&(rows.model=='random_forest')&(rows.subset=='all')].iloc[0]
        self.assertAlmostEqual(v.balanced_accuracy,balanced_accuracy_score(truth,p>=.5),places=12)

    def test_reliability_bins_account_for_every_observation(self):
        c=pd.read_csv(self.ext/'classification_reliability_seed0.csv');p=pd.read_csv(self.ext/'classification_test_predictions_seed0.csv')
        self.assertEqual(c.n.sum(),len(p))
        self.assertAlmostEqual((c.n*c.observed_frequency).sum(),p.y.sum(),places=8)
        self.assertAlmostEqual((c.n*c.mean_probability).sum(),p.probability.sum(),places=8)

    def test_sequential_model_independent_integrator(self):
        for k1,k2 in [(.015,.002),(.002,.015),(.01,.01),(.01,.0100000000001),(0.,.01),(.015,0.)]:
            matrix=np.array([[-k1,0,0],[k1,-k2,0],[0,k2,0]])
            for t in [0.,30.,120.,600.]:
                # DOP853 avoids a near-degenerate eigensystem reference problem.
                expected=solve_ivp(lambda _,y:matrix@y,(0,t),[1.,0.,0.],method='DOP853',rtol=1e-12,atol=1e-13).y[:,-1]
                actual=sequential_solution(t,k1,k2)
                np.testing.assert_allclose(actual,expected,atol=1e-10,rtol=1e-10)
                self.assertGreaterEqual(actual.min(),-1e-12);self.assertAlmostEqual(actual.sum(),1,places=12)

    def test_flow_current_includes_both_steps(self):
        # A->P->D: final P costs 2e, final D costs 4e.
        for q in [.2,.5,1.,2.]:
            row=flow_metrics(q);feed=q*1e-3*.1/60
            required=96485.33212*feed*(2*row['product_yield']+4*row['degraded_fraction'])
            self.assertAlmostEqual(required,row['minimum_current_A'],places=12)
            self.assertEqual(row['charge_feasible'],required<=.2)
        self.assertFalse(flow_metrics(1.)['charge_feasible'])
        self.assertGreater(flow_metrics(2.)['candidate_FE_at_available_current'],1)

    def test_flow_units_and_grid_constraints(self):
        grid=pd.read_csv(self.ext/'flow_grid.csv')
        for q in [.2,.5]:
            row=flow_metrics(q);moles_hour=q*1e-3*60*.1*row['product_yield']
            self.assertAlmostEqual(row['product_g_h'],moles_hour*200,places=12)
            self.assertAlmostEqual(row['sty_g_L_h'],row['product_g_h']/.001,places=10)
            self.assertAlmostEqual(row['energy_kWh_kg_at_assumed_3V'],.0006/(row['product_g_h']/1000),places=12)
        eligible=grid[grid.charge_feasible&(grid.conversion>=.8)&(grid.selectivity>=.8)]
        best=eligible.loc[eligible.sty_g_L_h.idxmax()]
        self.assertAlmostEqual(best.q_ml_min,.55,places=12)

    def test_molecular_36_jobs_and_6_optimizations(self):
        s=json.loads((self.mol/'summary.json').read_text());self.assertEqual(s['completed'],6);self.assertFalse(s['failures'])
        files=list(self.mol.glob('*/*/run.json'));self.assertEqual(len(files),36)
        for f in files:
            r=json.loads(f.read_text());self.assertTrue(r['completed']);self.assertEqual(r['returncode'],0)
            if f.parent.name=='neutral':self.assertTrue(r['optimization_converged'])
            j=json.loads((f.parent/'xtbout.json').read_text())
            self.assertAlmostEqual(j['total energy'],r['energy_Eh'],places=7)

    def test_molecular_fixed_nuclei_atom_identity_and_spin(self):
        for identity in self.mol.glob('*/identity.json'):
            d=json.loads(identity.read_text());base=Chem.AddHs(Chem.MolFromSmiles(d['input_smiles']))
            self.assertEqual(base.GetNumAtoms(),d['atom_count'])
            coord=(identity.parent/'neutral/xtbopt.xyz').read_bytes()
            symbols=[line.split()[0] for line in coord.decode().splitlines()[2:] if line.strip()]
            self.assertEqual(symbols,[a.GetSymbol() for a in base.GetAtoms()])
            for name,charge in [('cation',1),('anion',-1),('gas_neutral',0),('gas_cation',1),('gas_anion',-1)]:
                folder=identity.parent/name
                self.assertEqual((folder/'input.xyz').read_bytes(),coord)
                n_electrons=sum(a.GetAtomicNum() for a in base.GetAtoms())-charge
                run=json.loads((folder/'run.json').read_text());uhf=int(run['args'][run['args'].index('--uhf')+1])
                self.assertEqual((n_electrons-uhf)%2,0)

    def test_molecular_energy_and_charge_conservation(self):
        df=pd.read_csv(self.mol/'molecular_summary.csv');atoms=pd.read_csv(self.mol/'atom_charge_responses.csv')
        for row in df.itertuples():
            root=self.mol/row.molecule;energies={n:json.loads((root/n/'run.json').read_text())['energy_Eh'] for n in ['neutral','cation','anion']}
            self.assertAlmostEqual(row.fixed_nuclei_removal_gap_eV,(energies['cation']-energies['neutral'])*27.211386245988,places=10)
            a=atoms[atoms.molecule==row.molecule]
            for folder,charge,col in [('neutral',0,'q_neutral'),('cation',1,'q_cation'),('anion',-1,'q_anion')]:
                raw=np.loadtxt(root/folder/'charges');self.assertLess(abs(raw.sum()-charge),1e-6)
                np.testing.assert_allclose(raw,a[col],atol=1e-12)
            self.assertLess(abs(a.electron_removal_response.sum()-1),1e-6)
            self.assertLess(abs(a.electron_addition_response.sum()-1),1e-6)

if __name__=='__main__':unittest.main()
