import os
import unittest
from programs import PAML
import pandas as pd


class TestPAML(unittest.TestCase):
    input_path = f"/data/test/codon_aligned_seq_data.fas"
    output_path = "/data/test/paml.out"
    aux_dir = "/data/test/paml_aux/"
    tree_path = "/data/test/paml_tree.nwk"
    control_filepath = "/data/test/paml.ctl"

    def test_creation(self):
        prog = PAML()

    def test_exec(self):
        prog = PAML()
        prog.exec(
            input_path=self.input_path,
            output_path=self.output_path,
            aux_dir=self.aux_dir,
            additional_params={"NSsites": "3"},
            input_tree_path=self.tree_path,
            control_file_path=self.control_filepath
        )
        self.assertTrue(os.path.exists(self.output_path))

    def test_parse_output(self):
        prog = PAML()
        prog.exec(
            input_path=self.input_path,
            output_path=self.output_path,
            aux_dir=self.aux_dir,
            additional_params={"NSsites": "3"},
            input_tree_path=self.tree_path,
            control_file_path=self.control_filepath
        )
        result = prog.parse_output(self.output_path)
        NEB_positive_selection_analysis = pd.DataFrame.from_dict(result["NEB_positive_selection_analysis"])
        self.assertTrue(NEB_positive_selection_analysis.loc[NEB_positive_selection_analysis["position"] == 50, "is_significant"].values[0])
        self.assertTrue(abs(NEB_positive_selection_analysis.loc[NEB_positive_selection_analysis["position"] == 50, "mean(w)"].values[0]-4.705) < 0.1)
        self.assertTrue(not NEB_positive_selection_analysis.loc[NEB_positive_selection_analysis["position"] == 15, "is_significant"].values[0])
        self.assertTrue(abs(NEB_positive_selection_analysis.loc[NEB_positive_selection_analysis["position"] == 15, "p(w>1)"].values[0] - 0.913) < 0.1)
        self.assertTrue(result["ws_inference"][1]["prop"] <= 1)
        self.assertTrue(result["ws_inference"][1]["w"] < 0.4)
        self.assertTrue(int(result["duration(minutes)"]) < 2)

    def tearDown(self):
        if os.path.exists(self.output_path):
            os.remove(self.output_path)

