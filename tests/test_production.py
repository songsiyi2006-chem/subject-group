"""Checks for the second supplied pipeline, its repairs and added scientific audits."""
from pathlib import Path
import hashlib,importlib.util,json,re,unittest
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem,rdFreeSASA
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'production';A=P/'results/audit';R=P/'results/repaired'
spec=importlib.util.spec_from_file_location('production_audit',P/'scripts/audit_pipeline.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

class ProductionEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((R/'production_benchmark_results.json').read_text())
        cls.summary=json.loads((A/'audit_summary.json').read_text())
    def test_supplied_source_identity_and_failure_preserved(self):
        meta=json.loads((P/'source/source_record.json').read_text());raw=(P/'source/run_production_pipeline_original.py').read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(),meta['extracted_script_sha256'])
        execution=json.loads((P/'results/original_attempt/execution.json').read_text())
        self.assertEqual(execution['returncode'],1);self.assertFalse(execution['result_written'])
        self.assertIn("has no attribute 'calcSASA'",(P/'results/original_attempt/stderr.log').read_text())
    def test_repaired_modules_and_patch_are_explicit(self):
        self.assertTrue(all(t['status']=='completed' for t in self.data.values()))
        patch=(P/'source/repair.patch').read_text();self.assertIn('+    total_sasa = rdFreeSASA.CalcSASA',patch)
        self.assertIn('GetRvdw',patch);self.assertIn('GetDoubleProp("SASA")',patch)
    def test_source_bo_grid_and_observation_counts(self):
        space=pd.read_csv(R/'task_a_candidates.csv');obs=pd.read_csv(R/'task_a_observations.csv')
        self.assertEqual(len(space),288);self.assertEqual(len(obs),14)
        self.assertEqual(len(obs.drop_duplicates(subset=audit.FEATURES)),14)
        self.assertNotIn(14,space.current_mA_cm2.unique());self.assertNotIn(303.15,space.temp_K.unique())
        self.assertEqual(round(obs.observed_yield.max(),2),self.data['Task_A_Physical_BO']['optimal_yield'])
    def test_bo_formula_ties_independently(self):
        pool=pd.read_csv(R/'task_a_candidates.csv');ys=audit.latent_yield(pool)
        self.assertAlmostEqual(ys.max(),91.3865,places=10);self.assertEqual(int(np.isclose(ys,ys.max()).sum()),6)
        sub=pool[(pool.solvent=='MeCN')&(pool.current_mA_cm2==16)&(pool.temp_K==298.15)&(pool.electrolyte!='nBu4NOAc')]
        np.testing.assert_allclose(audit.latent_yield(sub),91.3865,atol=1e-10)
    def test_bo_paired_design_regret_and_no_duplicates(self):
        obs=pd.read_csv(A/'task_a_observations.csv');runs=pd.read_csv(A/'task_a_paired_runs.csv')
        self.assertEqual(len(obs),30*4*14);self.assertEqual(len(runs),120)
        for seed,g in obs.groupby('seed'):
            starts=[h.head(6)[['pool_index','observed']].to_numpy() for _,h in g.groupby('method')]
            for x in starts[1:]:np.testing.assert_array_equal(starts[0],x)
        for row in runs.itertuples():
            g=obs[(obs.seed==row.seed)&(obs.method==row.method)];self.assertEqual(g.pool_index.nunique(),14)
            best=g.loc[g.observed.idxmax()];self.assertAlmostEqual(91.3865-best.latent,row.recommendation_regret,places=10)
    def test_target_core_mapping_and_blocked_sites(self):
        mol=audit.make_molecule();self.assertEqual(mol.GetNumHeavyAtoms(),17);self.assertEqual(mol.GetNumAtoms(),30)
        self.assertEqual(mol.GetAtomWithIdx(6).GetSymbol(),'N')
        self.assertEqual({x.GetIdx() for x in mol.GetAtomWithIdx(7).GetNeighbors()},{6,8,11})
        self.assertEqual({x.GetIdx() for x in mol.GetAtomWithIdx(2).GetNeighbors()},{1,3,10})
        self.assertTrue(mol.GetAtomWithIdx(8).GetIsAromatic())
        for idx in [2,7]:self.assertFalse(any(n.GetAtomicNum()==1 for n in mol.GetAtomWithIdx(idx).GetNeighbors()))
        self.assertTrue(any(n.GetAtomicNum()==1 for n in mol.GetAtomWithIdx(8).GetNeighbors()))
        self.assertNotIn(15,audit.CORE_LABELS)
    def test_sasa_radii_and_single_geometry_partition(self):
        atoms=pd.read_csv(R/'task_b_all_atoms.csv');self.assertTrue((atoms.radius_A>0).all())
        self.assertAlmostEqual(atoms.sasa_A2.sum(),self.data['Task_B_3D_Cheminformatics']['total_sasa_A2'],places=9)
        self.assertEqual(self.summary['B']['classified_radii_unique'],[0.0])
        mol=Chem.MolFromMolFile(str(R/'task_b_mmff_geometry.mol'),removeHs=False)
        radii=[Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()]
        recalc=rdFreeSASA.CalcSASA(mol,radii)
        # Molfile coordinates have only four decimals: use an absolute area tolerance.
        self.assertAlmostEqual(recalc,atoms.sasa_A2.sum(),delta=.05)
    def test_conformer_convergence_partition_and_charge_invariance(self):
        conf=pd.read_csv(A/'task_b_conformers.csv');atoms=pd.read_csv(A/'task_b_conformer_atoms.csv')
        self.assertEqual(len(conf),32);self.assertTrue((conf.mmff_status==0).all())
        np.testing.assert_allclose(conf.total_sasa_A2,conf.atomic_sasa_sum,atol=1e-9)
        self.assertTrue((atoms.groupby('atom_index').gasteiger_charge.std()<1e-14).all())
        self.assertEqual(len(atoms),32*30)
    def test_conformer_site_scope_and_rank_count(self):
        sites=pd.read_csv(A/'task_b_conformer_sites.csv');summ=pd.read_csv(A/'task_b_site_summary.csv')
        self.assertEqual(sites.atom_index.nunique(),10);self.assertEqual(summ.atom_index.nunique(),9)
        self.assertEqual(int(summ.top_count.sum()),32);self.assertNotIn(0,summ.atom_index.to_numpy())
        self.assertEqual(int(summ[summ.site_label=='C3'].top_count.iloc[0]),0)
    def test_group_holdouts_do_not_leak_group(self):
        source=pd.read_csv(R/'task_c_all_sites.csv');pred=pd.read_csv(A/'task_c_predictions.csv')
        self.assertEqual(len(source),16)
        for row in pred.itertuples():self.assertEqual(str(source.loc[row.row_id,row.split]),str(row.held_out))
        metrics=pd.read_csv(A/'task_c_group_holdouts.csv')
        self.assertTrue((metrics.n_train==12).all());self.assertTrue((metrics.n_test==4).all())
        for row in metrics.itertuples():
            sub=pred[(pred.split==row.split)&(pred.held_out==row.held_out)&(pred.model==row.model)]
            self.assertAlmostEqual(np.mean(abs(sub.target-sub.prediction)),row.test_mae,places=10)
    def test_flow_source_agrees_with_independent_analytic_solution(self):
        original=self.data['Task_D_Flow_Electrochemistry']['flow_sweep_kinetics']
        for q in [100,300,600,1200]:
            row,profile=audit.flow_solution(q);o=original[f'{q}_uL_min']
            self.assertAlmostEqual(100*row['conversion'],o['conversion_pct'],delta=.0051)
            self.assertAlmostEqual(row['STY_mmol_L_h'],o['space_time_yield_mmol_L_h'],delta=.0051)
            self.assertLess(row['analytic_ode_error'],1e-8)
            self.assertTrue((profile.Cs_mol_m3>=0).all());self.assertTrue((profile.Cs_mol_m3<=profile.Cb_mol_m3).all())
            self.assertTrue((np.diff(profile.Cb_mol_m3)<=0).all())
    def test_flow_current_flux_closure(self):
        F=96485.33;Rgas=8.314;T=298.15
        for q in [100,300,600,1200]:
            row,p=audit.flow_solution(q)
            bv=.05*(p.Cs_mol_m3/50*np.exp(.5*F*.45/(Rgas*T))-np.exp(-.5*F*.45/(Rgas*T)))
            np.testing.assert_allclose(bv,p.j_A_m2,rtol=1e-12,atol=1e-12)
            lam=row['Da']/.1;ceq=row['Ceq_mol_m3'];keff=row['keff_m_s']
            current=quad(lambda z:F*keff*(50-ceq)*np.exp(-lam*z)*.01,0,.1,epsabs=1e-12)[0]
            self.assertAlmostEqual(current,row['current_from_molar_balance_A'],places=12)
    def test_flow_limiting_regimes_and_trend(self):
        zero,_=audit.flow_solution(300,eta=0);self.assertLess(abs(zero['conversion']),1e-12)
        slow,_=audit.flow_solution(300,j0=1e-10);fast,_=audit.flow_solution(300,j0=1e5)
        self.assertLess(slow['conversion'],1e-5);self.assertLess(abs(fast['keff_m_s']/fast['km_m_s']-1),1e-6)
        rows=[audit.flow_solution(q)[0] for q in [100,300,600,1200]]
        self.assertTrue(np.all(np.diff([r['km_m_s'] for r in rows])>0))
        self.assertTrue(np.all(np.diff([r['mass_transfer_resistance_fraction'] for r in rows])<0))
        self.assertTrue(np.all(np.diff([r['conversion'] for r in rows])<0))
    def test_target_xtb_records_and_charge_sums(self):
        s=self.summary['target_xTB'];self.assertEqual(s['status'],'completed')
        base=P/'results/target_xtb';geom=(base/'neutral/xtbopt.xyz').read_bytes()
        for folder,target in [('neutral',0),('cation',1),('anion',-1)]:
            run=json.loads((base/folder/'run.json').read_text());self.assertTrue(run['completed'])
            self.assertLess(abs(np.loadtxt(base/folder/'charges').sum()-target),1e-6)
            if folder!='neutral':self.assertEqual((base/folder/'input.xyz').read_bytes(),geom)
        atoms=pd.read_csv(base/'charge_responses.csv')
        self.assertAlmostEqual(atoms.removal_response.sum(),1,places=6)
        self.assertAlmostEqual(atoms.addition_response.sum(),1,places=6)
        self.assertAlmostEqual((s['energies_Eh']['cation']-s['energies_Eh']['neutral'])*27.211386245988,s['removal_difference_eV'],places=10)

if __name__=='__main__':unittest.main()
