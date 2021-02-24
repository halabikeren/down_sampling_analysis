import sys
import unittest
from tests import TestPDA, TestCDHIT, TestRate4ite, TestRandom, TestPAML, TestPipeline, TestSimulationPipeline

TESTS = [TestCDHIT, TestPDA, TestRandom, TestRate4ite, TestPAML, TestPipeline, TestSimulationPipeline]

suite = unittest.TestSuite()
loader = unittest.TestLoader()
for test in TESTS:
    suite.addTests(loader.loadTestsFromTestCase(test))
runner = unittest.TextTestRunner()
runner.run(suite)
ret = not runner.run(suite).wasSuccessful()
sys.exit(ret)
