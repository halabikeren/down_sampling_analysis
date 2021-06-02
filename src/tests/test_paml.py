import os
import unittest
from programs import Paml
import pandas as pd


class TestPAML(unittest.TestCase):
    input_path = f"/data/test/codon_aligned_seq_data.fas"
    output_path = "/data/test/paml.out"
    aux_dir = "/data/test/paml_aux/"

    def test_creation(self):
        prog = Paml()

    def test_exec(self):
        prog = Paml()
        prog.exec(
            input_path=self.input_path,
            output_path=self.output_path,
            aux_dir=self.aux_dir,
            additional_params={"NSsites": "3"},
        )
        self.assertTrue(os.path.exists(self.output_path))

    def test_parse_output(self):
        prog = Paml()
        prog.exec(
            input_path=self.input_path,
            output_path=self.output_path,
            aux_dir=self.aux_dir,
            additional_params={"NSsites": "3"},
        )
        result = prog.parse_output(self.output_path)
        neb_positive_selection_analysis = pd.DataFrame.from_dict(
            result["NEB_positive_selection_analysis"]
        )
        self.assertTrue(
            neb_positive_selection_analysis.loc[
                neb_positive_selection_analysis["position"] == 50, "is_significant"
            ].values[0]
        )
        self.assertTrue(
            abs(
                neb_positive_selection_analysis.loc[
                    neb_positive_selection_analysis["position"] == 50, "mean(w)"
                ].values[0]
                - 4.705
            )
            < 0.1
        )
        self.assertTrue(
            not neb_positive_selection_analysis.loc[
                neb_positive_selection_analysis["position"] == 15, "is_significant"
            ].values[0]
        )
        self.assertTrue(
            abs(
                neb_positive_selection_analysis.loc[
                    neb_positive_selection_analysis["position"] == 15, "p(w>1)"
                ].values[0]
                - 0.913
            )
            < 0.1
        )
        self.assertTrue(result["selection_parameters"][1]["prop"] <= 1)

        self.assertTrue(
            result["selection_parameters"][0]["w"]
            <= result["selection_parameters"][1]["w"]
        )
        self.assertTrue(
            result["selection_parameters"][1]["w"]
            <= result["selection_parameters"][2]["w"]
        )
        self.assertTrue(int(result["duration(minutes)"]) < 2)

    def tearDown(self):
        if os.path.exists(self.output_path):
            os.remove(self.output_path)
