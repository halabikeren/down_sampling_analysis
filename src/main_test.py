import unittest
import sys
from tests import (
    TestPDA,
    TestCDHIT,
    TestRate4ite,
    TestRandom,
    TestPAML,
    TestPipeline,
    TestSimulationPipeline,
)

TESTS = [
    TestCDHIT,
    TestPDA,
    TestRandom,
    TestRate4ite,
    TestPAML,
    TestPipeline,
    TestSimulationPipeline,
]


suite = unittest.TestSuite()
loader = unittest.TestLoader()
for test in TESTS:
    suite.addTests(loader.loadTestsFromTestCase(test))
runner = unittest.TextTestRunner()
ret = not runner.run(suite).wasSuccessful()
sys.exit(ret)
