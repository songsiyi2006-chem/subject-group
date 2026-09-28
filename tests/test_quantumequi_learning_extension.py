"""Small mathematical/data tests; never rerun the training study or quantum jobs."""
import copy
import math
import unittest

import numpy as np
import torch

from quantumequi.scripts import learning_extension as le


class LearningExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)
        cls.rows = le.read_reference()
        cls.prep = le.fit_preprocessing(cls.rows)

    def model(self, seed=11):
        torch.manual_seed(seed)
        return le.HydrogenDistanceMPNN(self.prep)

    def test_frozen_partition_exact_counts(self):
        self.assertEqual({s:sum(r["split"]==s for r in self.rows) for s in ["train","validation","test","ood_stretch"]},
                         {"train":17,"validation":8,"test":8,"ood_stretch":9})
        self.assertEqual(len({r["point_id"] for r in self.rows}), 42)

    def test_preprocessing_has_no_heldout_leakage(self):
        changed = copy.deepcopy(self.rows)
        for r in changed:
            if r["split"] != "train":
                r["energy_total_Hartree"] = 1e8
                r["dE_dR_Hartree_A"] = -1e9
                r["R_A"] = 500
        self.assertEqual(le.fit_preprocessing(changed), self.prep)

    def test_training_scaling_matches_independent_statistics(self):
        training = [r for r in self.rows if r["split"] == "train"]
        self.assertAlmostEqual(self.prep["distance_mean"], 1.14, places=14)
        self.assertAlmostEqual(self.prep["energy_std"], np.std([r["energy_total_Hartree"] for r in training]), places=15)

    def test_matched_initial_state_reproducible(self):
        a, b, c = self.model(40), self.model(40), self.model(41)
        self.assertEqual(le.state_digest(a.state_dict()), le.state_digest(b.state_dict()))
        self.assertNotEqual(le.state_digest(a.state_dict()), le.state_digest(c.state_dict()))

    def test_scalar_energy_shape_and_force_units_identity(self):
        model = self.model()
        e, g = le.radial_predictions(model, [0.7, 1.4])
        self.assertEqual(tuple(e.shape), (2,))
        self.assertEqual(tuple(g.shape), (2,))
        h = 1e-5
        with torch.no_grad():
            en = model(le.coordinates_from_distances([0.7+h,0.7-h,1.4+h,1.4-h])).numpy()
        np.testing.assert_allclose(g.detach().numpy(), (en[::2]-en[1::2])/(2*h), atol=1e-10)

    def test_complete_cartesian_force_difference(self):
        model = self.model()
        x = le.coordinates_from_distances([0.91])
        _, f = le.energy_forces(model, x)
        h=1e-5
        for k in range(6):
            xp,xm=x.clone(),x.clone()
            xp.reshape(-1)[k]+=h;xm.reshape(-1)[k]-=h
            with torch.no_grad():
                derivative=float((model(xp)-model(xm))[0])/(2*h)
            self.assertAlmostEqual(float(f.reshape(-1)[k]),-derivative,places=9)

    def test_translation_rotation_reflection_permutation(self):
        model=self.model()
        x=le.coordinates_from_distances([0.74,1.8])
        e,f=le.energy_forces(model,x)
        a=.49
        q=np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,-1.]])
        # Rotate original z axis to a mixed direction as well as reflect.
        q=q @ np.array([[0,0,1],[1,0,0],[0,1,0.]])
        transformed=(x.numpy()@q.T+[4,-3,2])[:,::-1].copy()
        en,fn=le.energy_forces(model,torch.tensor(transformed,dtype=torch.float64))
        np.testing.assert_allclose(en.detach(),e.detach(),atol=1e-14)
        np.testing.assert_allclose(fn.detach(),(f.detach().numpy()@q.T)[:,::-1],atol=1e-14)

    def test_force_and_torque_balance(self):
        model=self.model()
        x=le.coordinates_from_distances([.6,1.2,2.4])
        _,f=le.energy_forces(model,x)
        np.testing.assert_allclose(f.detach().numpy().sum(axis=1),0,atol=1e-15)
        np.testing.assert_allclose(np.cross(x.numpy(),f.detach().numpy()).sum(axis=1),0,atol=1e-15)

    def test_gradient_loss_backpropagates_to_message_parameters(self):
        model=self.model()
        _,g=le.radial_predictions(model,[.7,1.3],create_graph=True)
        loss=((g-torch.tensor([.2,-.3]))**2).mean()
        loss.backward()
        gradient=model.messages[0][0].weight.grad
        self.assertIsNotNone(gradient)
        self.assertTrue(torch.isfinite(gradient).all())
        self.assertGreater(float(torch.linalg.vector_norm(gradient)),1e-8)

    def test_energy_training_can_update_real_parameters(self):
        model=self.model()
        before=le.state_digest(model.state_dict())
        optimizer=torch.optim.Adam(model.parameters(),lr=.001)
        optimizer.zero_grad()
        energy=model(le.coordinates_from_distances([.7,1.3]))
        ((energy-torch.tensor([-1.1,-.95]))**2).mean().backward()
        optimizer.step()
        self.assertNotEqual(before,le.state_digest(model.state_dict()))

    def test_checkpoint_clone_detached_from_future_updates(self):
        model=self.model()
        state=copy.deepcopy(model.state_dict())
        digest=le.state_digest(state)
        with torch.no_grad():
            model.hydrogen_embedding.add_(1)
        self.assertEqual(digest,le.state_digest(state))
        model.load_state_dict(state)
        self.assertEqual(digest,le.state_digest(model.state_dict()))

    def test_common_validation_score_uses_energy_and_gradient(self):
        e=self.prep["energy_std"]*.2
        g=self.prep["gradient_std"]*.7
        self.assertAlmostEqual(le.validation_score(e,g,self.prep),.9)

    def test_ensemble_sample_standard_deviation(self):
        rows=[]
        for objective in le.OBJECTIVES:
            for split in ["train","validation","test","ood_stretch"]:
                for seed,e in zip([1,2,3],[-1.,0.,1.]):
                    rows.append({"objective":objective,"point_id":split,"split":split,"seed":seed,"R_A":1.,
                                 "predicted_energy_Hartree":e,"predicted_gradient_Hartree_A":2*e,
                                 "reference_energy_Hartree":0.,"reference_gradient_Hartree_A":0.})
        points,summary=le.ensemble_tables(rows)
        self.assertEqual(len(points),8)
        self.assertEqual(len(summary),8)
        self.assertEqual(points[0]["ensemble_energy_std_Hartree"],1.)
        self.assertEqual(points[0]["ensemble_gradient_std_Hartree_A"],2.)
        self.assertEqual(points[0]["ensemble_energy_abs_error_Hartree"],0.)
        with self.assertRaises(ValueError):
            le.ensemble_tables(rows[:-1])

    def test_invalid_geometry_rejected(self):
        model=self.model()
        with self.assertRaises(ValueError):
            le.coordinates_from_distances([0])
        with self.assertRaises(ValueError):
            model(torch.zeros((1,3,3),dtype=torch.float64))
        with self.assertRaises(ValueError):
            model(torch.zeros((1,2,3),dtype=torch.float64))
        with self.assertRaises(ValueError):
            le.coordinates_from_distances([float("nan")])


if __name__ == "__main__":
    unittest.main()
