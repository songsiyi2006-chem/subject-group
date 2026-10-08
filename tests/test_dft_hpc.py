"""Synthetic temporary fixtures test bookkeeping, never chemical accuracy."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

from dft_hpc import workflow as w


class DftHpcTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='synthetic_dft_test_')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / 'campaign'
        self.campaign = w.prepare(self.root)

    def save_campaign(self):
        w.write_json(self.root / 'campaign.json', self.campaign)

    def result(self, planned_job, attempt='synthetic_1', **changes):
        job = planned_job
        folder = self.root / 'runs' / job['id'] / attempt
        folder.mkdir(parents=True)
        raw = folder / 'engine.out'
        raw.write_text('SYNTHETIC UNIT TEST FIXTURE; NOT A PSI4 CALCULATION\n')
        record = dict(job=copy.deepcopy(job), code_sha256=self.campaign['code_sha256'],
                      status='completed', scf_converged=True,
                      energy_hartree=-10.0 if job['charge'] == 1 else -10.1,
                      engine_version='SYNTHETIC-TEST-ONLY',
                      geometry_sha256=job['geometry_sha256'],
                      engine_output_sha256=w.digest(raw), fixture_type='synthetic')
        record.update(changes)
        # Explicitly allow NaN here to test rejection of invalid external input.
        (folder / 'result.json').write_text(json.dumps(record), encoding='utf-8')
        return folder

    def collect(self):
        return w.collect(self.root / 'campaign.json')

    def test_prepare_four_jobs_and_refuse_overwrite(self):
        self.assertEqual(len(self.campaign['jobs']), 4)
        self.assertEqual(self.campaign['status'], 'prepared_not_executed')
        self.assertFalse((self.root / 'runs').exists())
        self.assertIn('--array=0-3%2', (self.root / 'submit.slurm').read_text())
        before = w.digest(self.root / 'campaign.json')
        with self.assertRaises(FileExistsError):
            w.prepare(self.root)
        self.assertEqual(before, w.digest(self.root / 'campaign.json'))

    def test_resources_rejected_before_creation(self):
        for kwargs in ({'cores': 0}, {'cores': True}, {'hours': 169},
                       {'memory_mb': 100}, {'concurrency': 17}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                w.prepare(Path(self.tmp.name) / 'invalid', **kwargs)
        self.assertFalse((Path(self.tmp.name) / 'invalid').exists())

    def test_xyz_and_electron_parity(self):
        p = Path(self.tmp.name) / 'synthetic.xyz'
        p.write_text('2\nsynthetic H2\nH 0 0 0\nH 0 0 0.74\n')
        atoms = w.read_xyz(p)
        self.assertEqual(w.validate_state(atoms, 0, 1), 2)
        self.assertEqual(w.validate_state(atoms, 1, 2), 1)
        for charge, multiplicity in [(0, 2), (1, 1), (3, 1), (0, 0), (True, 1)]:
            with self.subTest(state=(charge, multiplicity)), self.assertRaises(ValueError):
                w.validate_state(atoms, charge, multiplicity)
        for text in ['2\nx\nH 0 0 0\n', '1\nx\nH nan 0 0\n',
                     '2\nx\nH 0 0 0\nH 0 0 0.1\n', '1\nx\nXx 0 0 0\n']:
            p.write_text(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                w.read_xyz(p)

    def test_missing_results_are_not_predictions(self):
        report = self.collect()
        self.assertEqual(report['accepted_chemical_predictions'], 0)
        self.assertEqual(len(report['pairs']), 2)
        self.assertTrue(all(p['status'] == 'incomplete' for p in report['pairs']))
        self.assertTrue(all(p['energy_difference_eV'] is None for p in report['pairs']))

    def test_synthetic_pair_is_only_diagnostic(self):
        for job in self.campaign['jobs']:
            self.result(job)
        report = self.collect()
        self.assertEqual(report['accepted_chemical_predictions'], 0)
        for pair in report['pairs']:
            self.assertEqual(pair['status'], 'diagnostic_pair')
            self.assertAlmostEqual(pair['energy_difference_hartree'], -0.1)
            self.assertAlmostEqual(pair['energy_difference_eV'], -0.1 * w.HARTREE_TO_EV)
            self.assertIn('not validated', pair['interpretation'])

    def test_different_protocols_never_pair(self):
        for field, value in [('method', 'b3lyp'), ('basis', 'def2-svp'),
                             ('geometry_sha256', 'f' * 64), ('environment', 'water')]:
            with self.subTest(field=field):
                original = copy.deepcopy(self.campaign['jobs'][1])
                altered = copy.deepcopy(original)
                altered[field] = value
                self.assertNotEqual(w.protocol(original), w.protocol(altered))
        self.campaign['jobs'][1]['method'] = 'b3lyp'
        self.campaign['jobs'][3]['geometry_sha256'] = 'f' * 64
        self.save_campaign()
        for job in self.campaign['jobs']:
            self.result(job)
        with self.assertRaises(ValueError):
            self.collect()

    def test_duplicate_success_is_ambiguous(self):
        for job in self.campaign['jobs'][:2]:
            self.result(job)
        self.result(self.campaign['jobs'][0], attempt='synthetic_duplicate')
        pair = self.collect()['pairs'][0]
        self.assertEqual(pair['status'], 'incomplete')
        self.assertIn('ambiguous', [s['status'] for s in pair['states']])

    def test_invalid_records_rejected(self):
        job = self.campaign['jobs'][0]
        cases = [({'energy_hartree': float('nan')}, 'missing_or_nonfinite_energy'),
                 ({'energy_hartree': float('inf')}, 'missing_or_nonfinite_energy'),
                 ({'energy_hartree': True}, 'missing_or_nonfinite_energy'),
                 ({'scf_converged': False}, 'not_scf_completed'),
                 ({'scf_converged': None}, 'not_scf_completed'),
                 ({'status': 'failed'}, 'not_scf_completed'),
                 ({'engine_version': ''}, 'missing_engine_version'),
                 ({'geometry_sha256': 'changed'}, 'geometry_mismatch'),
                 ({'code_sha256': {}}, 'runner_identity_mismatch'),
                 ({'job': {}}, 'job_identity_mismatch')]
        for i, (changes, reason) in enumerate(cases):
            self.result(job, attempt=f'synthetic_{i}', **changes)
        report = self.collect()
        for audit, (_, expected) in zip(report['records'], cases):
            self.assertIn(expected, audit['reasons'])
        self.assertTrue(all(p['status'] == 'incomplete' for p in report['pairs']))

    def test_raw_log_changed_or_missing_rejected(self):
        first = self.result(self.campaign['jobs'][0])
        second = self.result(self.campaign['jobs'][1])
        (first / 'engine.out').write_text('SYNTHETIC TAMPERED LOG')
        (second / 'engine.out').unlink()
        self.assertTrue(all('raw_log_missing_or_changed' in r['reasons']
                            for r in self.collect()['records']))

    def test_engine_versions_must_match(self):
        self.result(self.campaign['jobs'][0])
        self.result(self.campaign['jobs'][1], engine_version='SYNTHETIC-OTHER')
        self.assertEqual(self.collect()['pairs'][0]['status'], 'incomplete')

    def test_malformed_record_is_audited(self):
        folder = self.result(self.campaign['jobs'][0])
        (folder / 'result.json').write_text('{broken')
        self.assertEqual(self.collect()['records'][0]['reasons'], ['invalid_result_record'])

    def test_nonobject_json_record_is_audited(self):
        folder = self.result(self.campaign['jobs'][0])
        for text in ['[]', 'null', '123', '"synthetic"']:
            with self.subTest(text=text):
                (folder / 'result.json').write_text(text)
                self.assertEqual(self.collect()['records'][0]['reasons'], ['invalid_result_record'])

    def test_empty_xyz_has_validation_error(self):
        p = Path(self.tmp.name) / 'empty.xyz'
        p.write_text('')
        with self.assertRaises(ValueError):
            w.read_xyz(p)

    def test_prepared_files_tampering_rejected(self):
        for filename in ['runner.py', self.campaign['jobs'][0]['geometry']]:
            p = self.root / filename
            original = p.read_bytes()
            p.write_bytes(original + b'\nSYNTHETIC TAMPER\n')
            with self.subTest(filename=filename), self.assertRaises(ValueError):
                self.collect()
            p.write_bytes(original)


if __name__ == '__main__':
    unittest.main()
