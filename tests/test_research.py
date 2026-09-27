"""Independent source/physical bookkeeping checks for the research extension."""
import ast,hashlib,importlib.util,json,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'research';A=R/'results/audit'
spec=importlib.util.spec_from_file_location('research_audit',R/'scripts/audit_research.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
class ResearchEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads((R/'results/original/research_grade_results.json').read_text(encoding='utf-8'))
        cls.summary=json.loads((A/'audit_summary.json').read_text(encoding='utf-8'))
    def test_archived_source_and_execution(self):
        record=json.loads((R/'source/source_record.json').read_text());run=json.loads((R/'results/original/execution.json').read_text())
        digest=hashlib.sha256((R/'source/run_research_engine_original.py').read_bytes()).hexdigest()
        self.assertEqual(digest,record['script_sha256']);self.assertEqual(digest,run['script_sha256']);self.assertEqual(run['returncode'],0)
        self.assertTrue(all(v['status']=='completed' for v in self.source.values()))
        tree=ast.parse((R/'source/run_research_engine_original.py').read_text(encoding='utf-8'))
        self.assertFalse(any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='GaussianProcessRegressor' for n in ast.walk(tree)))
    def test_pareto_independent_comparison_and_duplicates(self):
        df=pd.read_csv(A/'conditions_source_rounded.csv');self.assertEqual(len(df),135)
        for row in df.itertuples():
            dominated=any(o.Yield_pct>=row.Yield_pct and o.FE_pct>=row.FE_pct and o.SEC_kWh_kg<=row.SEC_kWh_kg and (o.Yield_pct>row.Yield_pct or o.FE_pct>row.FE_pct or o.SEC_kWh_kg<row.SEC_kWh_kg) for o in df.itertuples())
            self.assertEqual(row.pareto,not dominated)
        self.assertEqual(int(df.pareto.sum()),12);self.assertEqual(len(df[df.pareto].drop_duplicates(audit.OBJECTIVES)),10)
    def test_source_optimum_and_sample_values(self):
        df=pd.read_csv(A/'conditions_source_rounded.csv');out=self.source['Module1_MultiObjective_BO']
        for row in out['pareto_solutions_sample']+[out['recommended_green_compromise']]:
            selected=df[(df.j_mA_cm2==row['j_mA_cm2'])&(df.catalyst_mol_pct==row['catalyst_mol_pct'])&(df.stirring_rpm==row['stirring_rpm'])].iloc[0]
            for key,value in row.items():self.assertAlmostEqual(selected[key],value,places=9)
        best=df[df.pareto].loc[df[df.pareto].Green_Score.idxmax()];self.assertEqual(int(best.condition_id),53)
    def test_sec_via_charge_time_and_mass(self):
        df=pd.read_csv(A/'conditions_unrounded.csv')
        for r in df.itertuples():
            # Produce 1 mmol in a 1 cm2 cell: derive time and energy independently.
            amount=.001;charge=amount*2*96485.33/(r.FE_pct/100);current=r.j_mA_cm2*.001
            seconds=charge/current;energy_kwh=current*r.U_cell_V*seconds/3600000;mass_kg=amount*.223
            self.assertAlmostEqual(energy_kwh/mass_kg,r.SEC_kWh_kg,places=12)
    def test_weight_cases_and_feasibility(self):
        w=pd.read_csv(A/'scalarization_weights.csv');grid=pd.read_csv(A/'conditions_unrounded.csv').set_index('condition_id')
        self.assertEqual(len(w),231);np.testing.assert_allclose(w[['w_yield','w_fe','w_sec']].sum(axis=1),1)
        self.assertTrue(grid.loc[w.condition_id,'pareto'].all())
        for r in self.summary['multiobjective']['constraint_scenarios']:
            ok=(grid.Yield_pct>=r['min_yield'])&(grid.FE_pct>=r['min_FE'])&(grid.SEC_kWh_kg<=r['max_SEC']);self.assertEqual(int(ok.sum()),r['feasible_count'])
    def test_gp_group_predictions_and_metrics(self):
        p=pd.read_csv(A/'gp_predictions.csv');m=pd.read_csv(A/'gp_group_metrics.csv')
        for r in p.itertuples():self.assertEqual(r.condition_id//27,r.fold)
        self.assertEqual(len(p),270)
        for r in m.itertuples():
            sub=p[(p.fold==r.fold)&(p.model==r.model)];self.assertEqual(len(sub),27);self.assertEqual(r.n_train,108)
            self.assertAlmostEqual(np.abs(sub['pred_'+r.objective]-sub['true_'+r.objective]).mean(),r.MAE,places=10)
    def test_molecular_identities_sites_and_charge_independence(self):
        desc=pd.read_csv(A/'scope_descriptors.csv');sites=pd.read_csv(A/'scope_ch_sites.csv');self.assertEqual(len(desc),8)
        for r in desc.itertuples():
            mol=Chem.AddHs(Chem.MolFromSmiles(r.input_smiles));a=mol.GetAtomWithIdx(r.primary_atom)
            self.assertEqual(a.GetAtomicNum(),6);self.assertTrue(any(n.GetAtomicNum()==1 for n in a.GetNeighbors()))
            self.assertEqual(a.GetIsAromatic(),r.primary_aromatic);self.assertEqual(r.no_geometry_charge_difference,0)
            self.assertEqual(int(sites[sites.molecule==r.molecule].sort_values('gasteiger_charge').iloc[0].atom_index),r.primary_atom)
            self.assertEqual(Chem.GetFormalCharge(mol),0)
        caffeine=desc[desc.molecule=='caffeine'].iloc[0];self.assertEqual(caffeine.primary_hybridization,'SP3')
    def test_source_yields_are_exact_seeded_normal_draws(self):
        rng=np.random.RandomState(42)
        for row in self.source['Module2_Substrate_Scope']['scope_summary_table']:
            value=row['estimated_E_ox_V_vs_SCE'];mean,sd=(88,3) if value<1.45 else ((72,4) if value<1.85 else (42,6))
            self.assertEqual(round(float(rng.normal(mean,sd)),1),row['predicted_electrochemical_yield_pct'])
        seeded=pd.read_csv(A/'scope_random_yield_seeds.csv');self.assertEqual(len(seeded),8000)
    def test_mmff_convergence_records(self):
        d=pd.read_csv(A/'scope_descriptors.csv');c=pd.read_csv(A/'scope_conformers.csv')
        self.assertTrue((d.source_mmff_status==0).all());self.assertEqual(len(c),64);self.assertTrue((c.status==0).all())
        for name in audit.IDS:self.assertEqual(sum(m is not None for m in Chem.SDMolSupplier(str(R/'results/molecules'/name/'conformers.sdf'),removeHs=False)),8)
    def test_all_xtb_jobs_and_charge_balance(self):
        x=pd.read_csv(A/'scope_xtb.csv');self.assertEqual(len(x),8)
        for r in x.itertuples():
            folder=R/'results/molecules'/r.molecule
            for state,target in [('neutral',0),('cation',1)]:
                run=json.loads((folder/state/'run.json').read_text());self.assertTrue(run['completed'])
                self.assertAlmostEqual(np.loadtxt(folder/state/'charges').sum(),target,places=6)
                if state=='neutral':self.assertTrue(run['optimization_converged'])
            self.assertEqual((folder/'neutral/xtbopt.xyz').read_bytes(),(folder/'cation/input.xyz').read_bytes())
            self.assertAlmostEqual((r.cation_Eh-r.neutral_Eh)*27.211386245988,r.removal_difference_eV,places=9)
    def test_transport_source_rounding_and_volume_integral(self):
        for row in self.summary['transport']['rows']:
            phi=row['phi'];value=3*quad(lambda x:x*np.sinh(phi*x)/np.sinh(phi),0,1,epsabs=1e-13)[0]
            self.assertAlmostEqual(value,row['eta'],places=11)
            original=self.source['Module3_POP_Pore_Transport']['diffusion_kinetic_evaluation'][row['regime']]['transport_sweep']
            v=next(r for r in original if r['pellet_radius_um']==row['radius_um'])
            self.assertAlmostEqual(row['phi'],v['thiele_modulus_phi'],delta=.00051)
            self.assertAlmostEqual(row['eta'],v['internal_effectiveness_factor_eta'],delta=.000051)
    def test_fvm_convergence_and_conservation(self):
        c=pd.read_csv(A/'transport_grid_convergence.csv')
        for _,sub in c.groupby(['regime','radius_um']):self.assertTrue(np.all(np.diff(sub.sort_values('n_cells').abs_error)<0))
        for phi in [.065,1.3,9.8]:
            eta,x,c=audit.sphere_fvm(phi,320);flux_eta=3*(2*320)*(1-c[-1])/phi**2
            self.assertAlmostEqual(eta,flux_eta,delta=1e-8);self.assertTrue((c>0).all());self.assertTrue((c<=1).all());self.assertTrue((np.diff(c)>0).all())
    def test_transport_limits_thresholds_and_source_false_claim(self):
        self.assertAlmostEqual(audit.eta_sphere(1e-6),1,places=12);self.assertAlmostEqual(audit.eta_sphere(1e4),.00029997,places=12)
        for r in self.summary['transport']['thresholds']:self.assertAlmostEqual(audit.eta_sphere(r['phi']),r['target_eta'],places=10)
        row=next(r for r in self.summary['transport']['rows'] if r['regime']=='Microporous_POP' and r['radius_um']==100)
        self.assertGreater(row['eta'],.4)
    def test_external_film_balance(self):
        e=pd.read_csv(A/'external_film_scenarios.csv');t=pd.read_csv(A/'transport_audit.csv')
        for r in e.itertuples():
            internal=t[(t.regime==r.regime)&(t.radius_um==r.radius_um)].iloc[0]
            self.assertAlmostEqual(3*r.Bi*(1-r.Cs_over_Cbulk),internal.phi**2*r.eta_overall,places=9)
            self.assertLessEqual(r.eta_overall,internal.eta)
if __name__=='__main__':unittest.main()
