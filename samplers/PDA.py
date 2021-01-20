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

    @staticmethod
    def get_farthest_leaves(node: Tree, node_to_height: t.Dict[str, float]) -> t.Tuple[
        t.List[str], float]:
        """
        :param node: node under which to two farthest leaves need to be detected
        :param node_to_height: a map of node to its height in the tree
        :return: (1) a list of the two farthest leaves in the tree whose root is node
                 (2) the distance between the two leaves
        """
        if node.is_leaf():
            node_to_height[node.name] = 0
            return [], 0

        children_dists = []
        leaves_distances_from_node = []
        for child in node.get_children():
            farthest_leaves, leaves_dist = PDA.get_farthest_leaves(child, node_to_height)
            [child_height, farthest_leaf_from_child] = node_to_height[child]
            children_dists.append((farthest_leaves, leaves_dist))
            leaves_distances_from_node.append([child_height + child.dist, farthest_leaf_from_child])

        max_leaves_dist_from_child = max(children_dists, key=itemgetter(1))
        first_farthest_leaf = max(leaves_distances_from_node, key=itemgetter(0))
        leaves_distances_from_node.remove(first_farthest_leaf)
        second_farthest_leaf = max(leaves_distances_from_node, key=itemgetter(0))
        node_to_height[node.name] = first_farthest_leaf[0]

        max_leaves_dist = max_leaves_dist_from_child[1]
        most_distant_leaves = max_leaves_dist_from_child[0]
        if max_leaves_dist < first_farthest_leaf[0] + second_farthest_leaf[0]:
            max_leaves_dist = first_farthest_leaf[0] + second_farthest_leaf[0]
            most_distant_leaves = [first_farthest_leaf[1], second_farthest_leaf[1]]

        return most_distant_leaves, max_leaves_dist

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
        assert (k > 1, "required sample size must be at least 2")
        assert (k < len(self.tree.get_leaf_names()), f"required sample size must be smaller than the dataset size {len(self.tree.get_leaf_names())}")
        while self.sample_size < k:
            self.do_pd_step(is_weighted)
        return self.sample_subtree.get_leaf_names()
