import unittest
from tests import TestPipeline, TestPDA, TestCDHIT, TestRate4ite, TestRandom, TestPAML

TESTS = [TestCDHIT, TestPDA, TestRandom, TestRate4ite, TestPAML, TestPipeline]

suite = unittest.TestSuite()
loader = unittest.TestLoader()
for test in TESTS:
    suite.addTests(loader.loadTestsFromTestCase(test))
runner = unittest.TextTestRunner()
runner.run(suite)
