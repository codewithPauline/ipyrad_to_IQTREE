import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / 'scripts/workflow.py'
spec = importlib.util.spec_from_file_location('workflow', CLI)
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)

class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='radseq test ')
        self.work = Path(self.temp.name)
        self.alignment = self.work / 'input with spaces.phy'
        self.alignment.write_text('4 4\na ACGT\nb ACGA\nc ACCT\nd TCGA\n')
    def tearDown(self):
        self.temp.cleanup()
    def cli(self, *args, env=None):
        return subprocess.run([sys.executable, str(CLI), *map(str, args)],
                              text=True, capture_output=True, env=env)
    def test_demo(self):
        self.assertEqual(workflow.phylip(ROOT / 'examples/demo.phy'), {'samples': 8, 'sites': 1200})
    def test_bad_alignments(self):
        for data in ['4 4\na ACGT\na ACGA\nc ACCT\nd TCGA\n',
                     '4 4\na ACG\nb ACGA\nc ACCT\nd TCGA\n',
                     '5 4\na ACGT\nb ACGA\nc ACCT\nd TCGA\n',
                     '4 4\na NNNN\nb ACGA\nc ACCT\nd TCGA\n',
                     '4 4\na ACGZ\nb ACGA\nc ACCT\nd TCGA\n']:
            with self.subTest(data=data):
                self.alignment.write_text(data)
                with self.assertRaises(ValueError):
                    workflow.phylip(self.alignment)
    def test_dry_run_does_not_create_output(self):
        prefix = self.work / 'results/demo'
        result = self.cli('infer', self.alignment, '--out', prefix, '--dry-run')
        self.assertEqual(result.returncode, 0, result.stderr)
        command = json.loads(result.stdout)['command']
        self.assertEqual(command[2], str(self.alignment))
        self.assertEqual(command[command.index('-T')+1], '1')
        self.assertFalse(prefix.parent.exists())
    def test_real_subprocess_argument_forwarding(self):
        stub = self.work / 'fake iqtree'
        captured = self.work / 'args.json'
        stub.write_text('#!' + sys.executable + '\nimport json,sys,os\nopen(os.environ["CAPTURE"], "w").write(json.dumps(sys.argv[1:]))\n')
        stub.chmod(0o755)
        result = self.cli('infer', self.alignment, '--out', self.work/'result prefix',
                          '--binary', stub, '--threads', '2', env={**os.environ, 'CAPTURE': str(captured)})
        self.assertEqual(result.returncode, 0, result.stderr)
        args = json.loads(captured.read_text())
        self.assertEqual(args[args.index('-s')+1], str(self.alignment))
        self.assertEqual(args[args.index('-T')+1], '2')
        self.assertNotIn('--redo', args)
    def test_external_failure_propagates(self):
        stub = self.work / 'fail'
        stub.write_text('#!/bin/sh\nexit 17\n'); stub.chmod(0o755)
        result = self.cli('infer', self.alignment, '--out', self.work/'failed', '--binary', stub)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('17', result.stderr)
    def test_missing_input(self):
        self.assertNotEqual(self.cli('check', self.work/'absent.phy').returncode, 0)
    def test_invalid_threads(self):
        result = self.cli('infer', self.alignment, '--out', self.work/'demo', '--threads', '0', '--dry-run')
        self.assertNotEqual(result.returncode, 0)
    def test_trim_protects_existing_files(self):
        result = self.cli('trim', self.alignment, '--out', self.work, '--adapter', 'ACGT', '--dry-run')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('empty or new', result.stderr)
    def test_placeholder_adapter(self):
        result = self.cli('trim', self.alignment, '--out', self.work/'trim', '--adapter', 'YOUR_ADAPTER', '--dry-run')
        self.assertNotEqual(result.returncode, 0)
    def test_historical_paths_rejected(self):
        result = self.cli('assemble', ROOT/'docs/legacy/ipyrad_step1.txt', '--dry-run')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('historical paths', result.stderr)
    def test_assembly_steps(self):
        params = self.work / 'params.txt'; params.write_text('example params for command test')
        result = self.cli('assemble', params, '--steps', '321', '--dry-run')
        self.assertNotEqual(result.returncode, 0)
        result = self.cli('assemble', params, '--steps', '123', '--threads', '4', '--dry-run')
        self.assertEqual(result.returncode, 0, result.stderr)
        command = json.loads(result.stdout)['command']
        self.assertEqual(command[-4:], ['-s', '123', '-c', '4'])
    def test_qc_command(self):
        result = self.cli('qc', self.alignment, '--out', self.work/'qc', '--dry-run')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['command'][0], 'fastqc')

if __name__ == '__main__':
    unittest.main()
