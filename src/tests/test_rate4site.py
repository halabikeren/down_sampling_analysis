import os
import unittest
from programs import Rate4Site
import pandas as pd


class TestRate4ite(unittest.TestCase):
    input_path = f"/data/test/aligned_seq_data.fas"
    output_path = "/data/test/r4s.out"
    aux_dir = "/data/test/r4s_aux/"

    def test_creation(self):
        prog = Rate4Site()

    def test_exec(self):
        prog = Rate4Site()
        prog.exec(
            input_path=self.input_path,
            output_path=self.output_path,
            aux_dir=self.aux_dir,
            additional_params={"-zj": ""},
        )
        self.assertTrue(os.path.exists(self.output_path))

    def test_parsing(self):
        prog = Rate4Site()
        prog.exec(
            input_path=self.input_path,
            output_path=self.output_path,
            aux_dir=self.aux_dir,
            additional_params={"-zj": ""},
        )
        result = prog.parse_output(
            output_path=self.output_path, job_output_dir=self.aux_dir
        )
        self.assertAlmostEqual(result["alpha"], 4.43237)
        self.assertTrue(abs(result["log_likelihood"] - (-156.342)) < 0.1)
        rates_by_position = pd.DataFrame.from_dict(result["rate_by_position"])
        self.assertEqual(
            rates_by_position.loc[
                rates_by_position["position"] == 4, "sequence"
            ].values[0],
            "G",
        )
        self.assertTrue(
            abs(
                rates_by_position.loc[
                    rates_by_position["position"] == 4, "rate"
                ].values[0]
                - 0.7202
            )
            < 0.1,
        )

    def tearDown(self):
        if os.path.exists(self.output_path):
            os.remove(self.output_path)


if __name__ == "__main__":
    unittest.main()
