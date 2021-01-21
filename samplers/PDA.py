from copy import deepcopy
from ete3 import Tree
import typing as t
from pydantic import BaseModel
from operator import itemgetter


class PDA(BaseModel):
    tree: Tree
    taxon_to_weight: t.Dict[str, float] = None
    norm_factor: float = 1
    sample_size: int = 0
    sample_subtree: Tree = Tree()
    pd_score: float = 0

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def get_farthest_leaves(node: Tree, node_to_height: t.Dict[str, t.List[float, str]]) -> t.Tuple[t.List[str], float]:
        """
        :param node: node under which to two farthest leaves need to be detected
        :param node_to_height: a map of node
        to its height in the tree, represented both by the the distance from its farthest leaf and the name of that
        respective leaf
        :return: (1) a list of the two farthest leaves in the tree whose root is node (2) the
        distance between the two leaves
        """
        if node.is_leaf():
            node_to_height[node.name] = [0, ""]
            return [], 0
        children_diameters = []
        node_diameter_components = [[0, ""], [0, ""]]
        for child in node.get_children():
            children_diameters.append(PDA.get_farthest_leaves(child, node_to_height))
            [child_height, child_farthest_leaf] = node_to_height[child.name]
            child_dist_from_node = child_height + child.dist
            if child_dist_from_node > node_diameter_components[0][0]:
                node_diameter_components[1] = node_diameter_components[0]
                node_diameter_components[0] = [child_dist_from_node, child_farthest_leaf]
            elif child_dist_from_node > node_diameter_components[1]:
                node_diameter_components[1] = child_farthest_leaf

        [max_child_diameter_leaves, max_child_diameter] = max(children_diameters, key=itemgetter(1))
        node_diameter = node_diameter_components[0][0] + node_diameter_components[1][0]
        node_diameter_leaves = [node_diameter_components[0][1], node_diameter_components[1][1]]
        if node_diameter > max_child_diameter:
            return node_diameter_leaves, node_diameter
        else:
            return max_child_diameter_leaves, max_child_diameter

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
        if self.sample_size == 0:
            farthest_leaves, pd_score = self.get_farthest_leaves(self.tree, dict())
            self.sample_subtree = deepcopy(self.tree)
            self.sample_subtree.prune(farthest_leaves)
            self.pd_score = self.norm_factor * pd_score
            if is_weighted:
                self.pd_score += sum([self.taxon_to_weight[taxon] for taxon in farthest_leaves])

        else:
            candidates = [node.name for node in self.tree.traverse() if
                          node.name not in [sub_node.name for sub_node in self.sample_subtree.traverse()]]
            max_pd_addition = 0
            max_candidate = max_pd_addition_attachment_node = None
            for candidate in candidates:
                candidate_parent_in_sample, pd_addition = self.get_pd_addition(candidate)
                if pd_addition > max_pd_addition:
                    max_pd_addition = pd_addition
                    max_candidate = candidate
                    max_pd_addition_attachment_node = candidate_parent_in_sample
            sample_parent = self.sample_subtree.search_nodes(name=max_pd_addition_attachment_node)[0]
            sample_parent.add_child(name=max_candidate, dist=max_pd_addition)
            self.pd_score += self.norm_factor * max_pd_addition
            if is_weighted:
                self.pd_score += self.taxon_to_weight[max_candidate]

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

    def compute_sample(self, k: int, is_weighted: bool = False) -> t.List[str]:
        """
        computes the most phylogenetically diverse weighted sample based on the greedy algorithm of Steel (2005).
        for more details see https://academic.oup.com/sysbio/article-abstract/54/4/527/2842877
        :param k: required sample size
        :param is_weighted: indicates weather the computed PD should be weighted or not
        :return: list of names of chosen leaves
        """
        if k == 1:
            raise ValueError("required sample size must be at least 2")
        elif k >= len(self.tree.get_leaf_names()):
            raise ValueError(
                f"required sample size must be smaller than the dataset size {len(self.tree.get_leaf_names())}")
        self.add_internal_names()
        while self.sample_size < k:
            self.do_pd_step(is_weighted)
        return self.sample_subtree.get_leaf_names()
