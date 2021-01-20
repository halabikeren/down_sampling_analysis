from ete3 import Tree
import typing as t
import pydantic

@pydantic
class PDA:
    tree: Tree
    taxon_to_weight: t.Dict[str, float] = None
    norm_factor: float = 1
    sample_size: int = 0
    sample_subtree: Tree = Tree()
    pd_score: float = 0

    @staticmethod
    def __init__(self, tree: Tree, taxon_to_weight: t.Dict[str, float] = None):
        self.tree = tree
        self.taxon_to_weight = taxon_to_weight

    def get_farthest_leaves(self) -> t.List[str]:
        """
        :return: returns a list with the names of the two farthest leaves in the tree
        """
        pass

    def get_pd_addition(self, leaf_name) -> t.Tuple[str, float, float]:
        """
        :return: (1) name of the parent of the leaf that is included in the current sample subtree
        (2) length of the branch that would be added to the sample subtree if the given leaf name is appended to the
        sample which equals to the addition to the non-weighted pd score as a result.
        """

    def do_pd_step(self, is_weighted: bool = False):
        """
        performs a PD step by addition of a leaf to the sample subtree and update of the pd score
        :param is_weighted: indicates weather the computed PD should be weighted or not
        :return: None
        """
        pass

    def compute_sample(self, is_weighted: bool = False):
        """
        computes the most phylogenetically diverse weighted sample based on the greedy algorithm of Steel (2005).
        for more details see https://academic.oup.com/sysbio/article-abstract/54/4/527/2842877
        :param is_weighted: indicates weather the computed PD should be weighted or not
        :return: None
        """
        pass
