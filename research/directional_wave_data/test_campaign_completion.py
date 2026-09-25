"""Prevent known false-completion modes; these tests never run OW3D or send mail."""
import pathlib
import tempfile
import unittest
from run_ow3d_campaign import native_completion

ROOT=pathlib.Path(__file__).resolve().parents[2]/'artifacts/campaign_controller_tests'
ROOT.mkdir(parents=True,exist_ok=True)

class CompletionTests(unittest.TestCase):
    def setUp(self):
        self.path=pathlib.Path(tempfile.mkdtemp(dir=ROOT))
        self.case={'expected_final_time_s':220}
        (self.path/'ow3d.log').write_text('JOB IS COMPLETE\n')
        (self.path/'OceanWave3D.end').write_text('header\n1 1 1025 257 2.200000000000004D+02\n')

    def test_successful_native_close(self):
        native_completion(self.case,self.path,0)

    def test_zero_exit_without_completion_is_failure(self):
        (self.path/'ow3d.log').write_text('Error reading input\n')
        with self.assertRaises(RuntimeError):native_completion(self.case,self.path,0)

    def test_wrong_final_time_is_failure(self):
        (self.path/'OceanWave3D.end').write_text('header\n1 1 1025 257 0.4\n')
        with self.assertRaises(RuntimeError):native_completion(self.case,self.path,0)

    def test_nonzero_exit_is_failure(self):
        with self.assertRaises(RuntimeError):native_completion(self.case,self.path,1)

if __name__=='__main__':unittest.main()
