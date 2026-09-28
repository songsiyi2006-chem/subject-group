"""Numerical/representation tests; none constitute chemistry validation."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFreeSASA

from electrograph.scripts.structure_reviewed import (
    MAPPED_EXAMPLE, R_KCAL_MOL_K, SOURCE_REACTION, aggregate_minima,
    audit_target_identity, heavy_atom_rmsd, mapped_bond_changes, minima_weights,
    parse_molecule, profile_molecule, radius_of_gyration, reaction_inventory,
    sasa_rotation_average, select_minima, vdw_radii,
)


class StructureReviewTests(unittest.TestCase):
    def test_weights_normalized_and_energy_shift_invariant(self):
        first = minima_weights([-50.0, -49.0, -45.0], 298.15)
        second = minima_weights([500.0, 501.0, 505.0], 298.15)
        np.testing.assert_allclose(first, second, atol=1e-15)
        self.assertAlmostEqual(sum(first), 1.0, places=14)
        self.assertGreater(first[0], first[1])

    def test_energy_units_analytic_two_state_ratio(self):
        w = minima_weights([0.0, R_KCAL_MOL_K * 300.0], 300.0)
        self.assertAlmostEqual(w[1]/w[0], np.exp(-1.0), places=14)
        self.assertAlmostEqual(R_KCAL_MOL_K*4184, 8.31446261815324, places=12)

    def test_weights_reject_invalid_energy_temperature(self):
        for e, t in [([], 300), ([np.nan], 300), ([[1]], 300), ([1], 0), ([1], np.inf)]:
            with self.assertRaises(ValueError):
                minima_weights(e, t)

    def test_high_temperature_flattens_populations(self):
        cool, hot = minima_weights([0, 1, 2], 250), minima_weights([0, 1, 2], 350)
        self.assertGreater(cool[0], hot[0])
        self.assertGreater(1/np.sum(hot**2), 1/np.sum(cool**2))

    def test_mass_weighted_rg_analytic(self):
        # Two masses, separation 2: Rg^2 = m1*m2*d^2/(m1+m2)^2 = 3/4.
        self.assertAlmostEqual(radius_of_gyration([[0, 0, 0], [2, 0, 0]], [1, 3]), np.sqrt(.75))
        self.assertAlmostEqual(radius_of_gyration([[0, 0, 0], [2, 0, 0]], [1, 1]), 1)

    def test_rg_translation_rotation_and_mass_scale_invariance(self):
        xyz = np.array([[0., 0, 0], [2, 1, -1], [-1, .5, 2]])
        rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
        self.assertAlmostEqual(radius_of_gyration(xyz, [1, 12, 16]),
            radius_of_gyration(xyz @ rotation + [25, -41, 5], [3, 36, 48]), places=13)
        with self.assertRaises(ValueError):
            radius_of_gyration(xyz, [1, 0, 1])

    def test_explicit_vdw_radii_and_sasa_bounded_by_isolated_spheres(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(mol, randomSeed=81)
        radii = vdw_radii(mol)
        self.assertGreater(min(radii), 0)
        options = rdFreeSASA.SASAOpts(rdFreeSASA.SASAAlgorithm.ShrakeRupley,
                                    rdFreeSASA.SASAClassifier.Protor, 1.4)
        area = rdFreeSASA.CalcSASA(mol, radii, opts=options)
        self.assertGreater(area, 0)
        self.assertLess(area, 4*np.pi*np.sum((np.array(radii)+1.4)**2))
        surface, values = sasa_rotation_average(mol, 0, radii)
        self.assertEqual(len(values), 24)
        self.assertGreaterEqual(surface["sasa_A2"], min(values))
        self.assertLessEqual(surface["sasa_A2"], max(values))
        self.assertGreaterEqual(surface["sasa_orientation_sd_A2"], 0)

    def test_rmsd_does_not_change_coordinates_and_dedup_prefers_lower_energy(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
        AllChem.EmbedMolecule(mol, randomSeed=8)
        duplicate = Chem.Conformer(mol.GetConformer())
        mol.AddConformer(duplicate, assignId=True)
        before = mol.GetConformer(0).GetPositions().copy()
        self.assertLess(heavy_atom_rmsd(mol, 0, 1), 1e-6)
        np.testing.assert_array_equal(before, mol.GetConformer(0).GetPositions())
        rows = [{"conformer_id": 0, "optimization_status": 0, "energy_kcal_mol": 1},
                {"conformer_id": 1, "optimization_status": 0, "energy_kcal_mol": 0}]
        selected, mapping = select_minima(mol, rows)
        self.assertEqual(selected, [1])
        self.assertTrue(all(r["representative_id"] == 1 for r in mapping))

    def test_unconverged_geometry_excluded_from_minima(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
        AllChem.EmbedMolecule(mol, randomSeed=9)
        selected, mapping = select_minima(mol, [{"conformer_id": 0, "optimization_status": 1, "energy_kcal_mol": -1000}])
        self.assertEqual(selected, [])
        self.assertEqual(mapping, [])
        self.assertEqual(aggregate_minima([], [], 298.15)[0]["status"], "no_converged_minima")

    def test_same_seed_reproduces_small_ensemble_and_coordinate_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = profile_molecule("CCO", "ethanol", "ethanol", Path(tmp)/"a", seeds=(41,), attempts=4)
            second = profile_molecule("CCO", "ethanol", "ethanol", Path(tmp)/"b", seeds=(41,), attempts=4)
            self.assertEqual(first["records"], second["records"])
            self.assertEqual(first["coordinate_sha256"], second["coordinate_sha256"])
            self.assertEqual(first["converged"]+first["nonconverged"], first["returned_after_embedding_pruning"])

    def test_unsupported_and_invalid_molecules_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(profile_molecule("[Fe]", "metal", "metal", Path(tmp), seeds=(1,), attempts=1)["status"], "unsupported")
        with self.assertRaises(ValueError):
            parse_molecule("")
        with self.assertRaises(ValueError):
            parse_molecule("C1(")

    def test_original_reaction_retains_formula_defect(self):
        record = reaction_inventory(SOURCE_REACTION)
        self.assertFalse(record["elementally_balanced"])
        self.assertFalse(record["heavy_atom_count_equal"])
        self.assertEqual(record["element_delta_products_minus_reactants"]["C"], -1)
        with self.assertRaises(ValueError):
            mapped_bond_changes(SOURCE_REACTION)

    def test_equal_heavy_counts_do_not_prove_elemental_balance(self):
        record = reaction_inventory("CC>>CO")
        self.assertTrue(record["heavy_atom_count_equal"])
        self.assertFalse(record["elementally_balanced"])
        self.assertEqual(record["element_delta_products_minus_reactants"], {"C": -1, "H": -2, "O": 1})

    def test_balanced_mapped_example_recovers_exact_four_bond_edits(self):
        record = mapped_bond_changes(MAPPED_EXAMPLE)
        self.assertTrue(record["inventory"]["elementally_balanced"])
        self.assertTrue(record["inventory"]["charge_balanced"])
        observed = {(tuple(c["atom_maps"]), c["kind"]) for c in record["bond_changes"]}
        self.assertEqual(observed, {((2, 3), "order_changed"), ((2, 4), "broken"),
                                    ((3, 5), "broken"), ((4, 5), "formed")})

    def test_reject_duplicate_missing_maps_and_bad_molecules(self):
        for rxn in ["[CH3:1][CH3:1]>>[CH3:1][CH3:2]", "CC>>CC",
                    "[CH3:1][CH3:2]>>[CH3:1][CH3:3]", "C1(>>CC", "CC>CC"]:
            with self.assertRaises(ValueError):
                mapped_bond_changes(rxn)

    def test_map_identity_mismatch_element_or_isotope_rejected(self):
        for rxn in ["[CH4:1]>>[OH2:1]", "[13CH4:1]>>[12CH4:1]"]:
            with self.assertRaises(ValueError):
                mapped_bond_changes(rxn)

    def test_implicit_hydrogen_and_charge_changes_are_audited(self):
        unbalanced = mapped_bond_changes("[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]")
        self.assertEqual(unbalanced["inventory"]["element_delta_products_minus_reactants"]["H"], -2)
        charged = mapped_bond_changes("[NH4+:1]>>[NH3:1]")
        self.assertEqual(charged["inventory"]["charge_delta"], -1)
        self.assertEqual(charged["atom_property_changes"][0]["atom_map"], 1)

    def test_target_numbering_confirms_five_methoxy_two_phenyl(self):
        identity = audit_target_identity()
        self.assertTrue(all(identity["checks"].values()))
        self.assertEqual(identity["formula"], "C15H13NO")
        self.assertEqual(identity["source_zero_based_atom_index_by_indole_position"]["5"], 2)


if __name__ == "__main__":
    unittest.main()
