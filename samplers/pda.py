import os
import re
import subprocess
import typing as t
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from operator import itemgetter

from Bio import SeqIO, AlignIO
from ete3 import Tree

from samplers import Sampler

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)


@dataclass
class Pda(Sampler):
    taxon_to_weight: t.Dict[str, float] = None
    taxon_to_weight_filepath: str = None
    norm_factor: float = 1
    sample_subtree: Tree = (
        Tree()
    )  # sample subtree is saved such that if later on, a larger sample is required, we can simply extend the current one rather than start over
    pd_score: float = 0

    def compute_taxon_weights(self, aligned_sequences_path: str):
        """
        computes the weight of each taxon (or sequence) based on the
        :return:
        """
        self.taxon_to_weight = dict()  # reset taxon to weight
        alignment = list(AlignIO.parse(aligned_sequences_path, "fasta"))[0]
        alignment_records = [record for record in alignment]
        # compute for each position the frequency of characters in it
        pos_to_char_feq = dict()
        for pos in range(alignment.get_alignment_length()):
            char_freq = defaultdict(float)
            for seq_record in alignment_records:
                char = seq_record.seq[pos]
                char_freq[char] += 1
            for char in char_freq:
                char_freq[char] /= len(alignment)
            pos_to_char_feq[pos] = char_freq
        # compute the weight of each taxon
        for seq_record in alignment_records:
            weight = 0
            for pos in pos_to_char_feq:
                pos_weight = 1 - pos_to_char_feq[pos]["-"]
                weight += pos_weight * pos_to_char_feq[pos][seq_record.seq[pos]]
            self.taxon_to_weight[
                seq_record.description if seq_record.description else seq_record.id
            ] = (weight / alignment.get_alignment_length())

    @staticmethod
    def get_max_pd_pair(
        node: Tree, node_to_height: t.Dict[str, t.Tuple[float, str]]
    ) -> t.Tuple[t.List[str], float]:
        """
        :param node: node under which to two farthest leaves need to be detected
        :param node_to_height: a map of node
        to its height in the tree, represented both by the the distance from its farthest leaf and the name of that
        respective leaf
        :return: (1) a list of the two farthest leaves in the tree whose root is node (2) the
        distance between the two leaves
        """
        if node.is_leaf():
            node_to_height[node.name] = (0, node.name)
            return [], 0
        children_diameters = []
        node_diameter_components = [(0, ""), (0, "")]
        for child in node.get_children():
            children_diameters.append(Pda.get_max_pd_pair(child, node_to_height))
            [child_height, child_farthest_leaf] = node_to_height[child.name]
            child_dist_from_node = child_height + child.dist
            if child_dist_from_node > node_diameter_components[0][0]:
                node_diameter_components[1] = node_diameter_components[0]
                node_diameter_components[0] = (
                    child_dist_from_node,
                    child_farthest_leaf,
                )
            elif child_dist_from_node > node_diameter_components[1][0]:
                node_diameter_components[1] = (
                    child_dist_from_node,
                    child_farthest_leaf,
                )

        node_to_height[node.name] = node_diameter_components[0]
        [max_child_diameter_leaves, max_child_diameter] = max(
            children_diameters, key=itemgetter(1)
        )
        node_diameter = node_diameter_components[0][0] + node_diameter_components[1][0]
        node_diameter_leaves = [
            node_diameter_components[0][1],
            node_diameter_components[1][1],
        ]
        if node_diameter > max_child_diameter:
            return node_diameter_leaves, node_diameter
        else:
            return max_child_diameter_leaves, max_child_diameter

    def get_max_wpd_pair(self) -> t.Tuple[t.List[str], float]:
        """
        :return: returns a pair of leaves maximizing the weighted pd score
        """
        leaves = self.tree.get_leaves()
        pairs = [
            (leaves[i], leaves[j])
            for i in range(len(leaves))
            for j in range(i + 1, len(leaves))
        ]
        max_pd_score = 0
        max_pd_leaves = ["", ""]
        for (leaf_1, leaf_2) in pairs:
            pd_score = self.norm_factor * self.tree.get_distance(leaf_1, leaf_2) + (
                self.taxon_to_weight[leaf_1.name] + self.taxon_to_weight[leaf_2.name]
            )
            if pd_score > max_pd_score:
                max_pd_score = pd_score
                max_pd_leaves = [leaf_1.name, leaf_2.name]
        return max_pd_leaves, max_pd_score

    def get_pd_addition(self, leaf_name) -> t.Tuple[str, float]:
        """
        :return: (1) name of the parent of the leaf that is included in the current sample subtree
        (2) length of the branch that would be added to the sample subtree if the given leaf name is appended to the
        sample which equals to the addition to the non-weighted pd score as a result.
        """
        branch_length = 0
        node = self.tree.search_nodes(name=leaf_name)[0]
        while len(self.sample_subtree.search_nodes(name=node.name)) == 0:
            branch_length += node.dist
            node = node.up

        return node.name, branch_length

    def do_pd_step(self, is_weighted: bool = False):
        """
        performs a PD step by addition of a leaf to the sample subtree and update of the pd score
        :param is_weighted: indicates weather the computed PD should be weighted or not
        :return: None
        """
        if len(self.sample_subtree.get_leaf_names()) < 2:
            if not is_weighted:
                chosen_leaves, pd_score = self.get_max_pd_pair(self.tree, dict())
            else:
                chosen_leaves, pd_score = self.get_max_wpd_pair()
            self.sample_subtree = deepcopy(self.tree)
            self.sample_subtree.prune(chosen_leaves)
            self.pd_score = self.norm_factor * pd_score
            if is_weighted:
                self.pd_score += sum(
                    [self.taxon_to_weight[taxon] for taxon in chosen_leaves]
                )
        else:
            candidates = [
                leaf
                for leaf in self.tree.get_leaf_names()
                if leaf not in self.sample_subtree.get_leaf_names()
            ]
            max_pd_addition = 0
            max_wpd_addition = 0
            max_candidate = max_pd_addition_attachment_node = None
            for candidate in candidates:
                candidate_parent_in_sample, pd_addition = self.get_pd_addition(
                    candidate
                )
                wpd_addition = self.norm_factor * pd_addition
                if is_weighted:
                    wpd_addition += self.taxon_to_weight[candidate]
                if wpd_addition > max_wpd_addition:
                    max_pd_addition = pd_addition
                    max_wpd_addition = wpd_addition
                    max_candidate = candidate
                    max_pd_addition_attachment_node = candidate_parent_in_sample
            sample_parent = self.sample_subtree.search_nodes(
                name=max_pd_addition_attachment_node
            )[0]
            sample_parent.add_child(name=max_candidate, dist=max_pd_addition)
            self.pd_score += max_wpd_addition

    def add_internal_names(self):
        """
        Adds names to internal nodes
        :return: None
        """
        i = 1
        for node in self.tree.traverse("postorder"):
            if node.name == "":
                node.name = "N" + str(i)
                i += 1

    def exec_external_pda(self, k: int, is_weighted: bool = False):
        if not os.path.exists(self.aux_dir):
            subprocess.run(f"mkdir -p {self.aux_dir}", shell=True, capture_output=True)
        self.tree.write(outfile=f"{self.aux_dir}/tree.nwk", format=5)
        weights_arg = ""
        if is_weighted:
            with open(
                f"{os.path.dirname(self.sequences_path)}/weights.txt", "w"
            ) as weights_file:
                weights_file.write(f"{self.norm_factor}\n")
                for taxon in self.taxon_to_weight:
                    weights_file.write(f"{taxon}\t{self.taxon_to_weight[taxon]}\n")
            weights_arg += f" -e {os.path.dirname(self.sequences_path)}/weights.txt"

        process = subprocess.run(
            f"{os.environ['pda']} -g -k {k}{weights_arg} {self.aux_dir}/tree.nwk {self.aux_dir}/out.pda",
            shell=True,
            capture_output=True,
        )
        if process.returncode != 0:
            raise RuntimeError(
                f"PDA failed to properly execute and provide an output file. Execution "
                f"output is {process.stderr}"
            )

        output_regex = re.compile(
            "For k = \d* the optimal PD score is (\d*).*?The optimal PD set has \d* taxa\:(.*?)Corresponding",
            re.MULTILINE | re.DOTALL,
        )
        with open(f"{self.aux_dir}/out.pda", "r") as result:
            result_content = result.read()
            output = output_regex.search(result_content)
        self.pd_score = float(output.group(1))
        sample_members = [
            member for member in output.group(2).split("\n") if member != ""
        ]
        self.sample_subtree = deepcopy(self.tree)
        self.sample_subtree.prune(sample_members)

        if self.aux_dir != os.path.dirname(os.path.realpath(self.sequences_path)):
            subprocess.run(f"rm -r {self.aux_dir}", shell=True, capture_output=True)

    def get_sample(
        self,
        k: int,
        is_weighted: bool = False,
        use_external=False,
    ) -> t.Union[str, t.List[SeqIO.SeqRecord]]:
        """
        computes the most phylogenetically diverse weighted sample based on the greedy algorithm of Steel (2005).
        for more details see https://academic.oup.com/sysbio/article-abstract/54/4/527/2842877
        :param k: required sample size
        :param is_weighted: indicates weather the computed PD should be weighted or not
        :param use_external: indicates weather the pda tool should be used or internally implemented code
        :return: list of names of chosen leaves
        """
        if k < 1:
            logger.error("Required sample size must be at least 2")
            raise ValueError("Required sample size must be at least 2")

        sample = super(Pda, self).get_sample(k)
        if k == len(self.sequences):
            self.sample_subtree = deepcopy(self.tree)
            self.pd_score = sum([node.dist for node in self.sample_subtree.traverse()])
            if is_weighted:
                self.pd_score = sum(
                    [
                        self.taxon_to_weight[taxon]
                        for taxon in self.sample_subtree.get_leaf_names()
                    ]
                )
            return sample

        else:
            if not use_external:
                self.add_internal_names()
                while len(self.sample_subtree.get_leaf_names()) < k:
                    self.do_pd_step(is_weighted=is_weighted)
            else:
                self.exec_external_pda(k, is_weighted=is_weighted)

            return [
                record
                for record in self.sequences
                if record.description in self.sample_subtree.get_leaf_names()
            ]
