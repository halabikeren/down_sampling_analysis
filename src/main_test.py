import unittest
from tests import TestPipeline, TestPDA, TestCDHIT, TestRate4ite, TestRandom

TESTS = (
    TestPipeline,
    TestPDA,
    TestCDHIT,
    TestRate4ite,
    TestRandom,
)

suite = unittest.TestSuite()
loader = unittest.TestLoader()
for test in TESTS:
    suite.addTests(loader.loadTestsFromTestCase(test))
runner = unittest.TextTestRunner()
runner.run(suite)
