"""Python wrappers to random tree generation in the IQ-TREE library."""

from enum import Enum, auto

from cogent3 import make_tree
from cogent3.core.tree import PhyloNode

from piqtree._piqtree import iq_random_tree
from piqtree.iqtree._decorator import iqtree_func
from piqtree.util import process_rand_seed_nonzero

iq_random_tree = iqtree_func(iq_random_tree)


class TreeGenMode(Enum):
    """Setting under which to generate a random tree."""

    YULE_HARDING = auto()
    UNIFORM = auto()
    CATERPILLAR = auto()
    BALANCED = auto()
    BIRTH_DEATH = auto()
    STAR_TREE = auto()


def _get_tree_gen_mode(tree_mode: TreeGenMode | str) -> TreeGenMode:
    if isinstance(tree_mode, TreeGenMode):
        return tree_mode
    try:
        return TreeGenMode[tree_mode.upper()]
    except KeyError:
        allowed = ", ".join(repr(mode.name) for mode in TreeGenMode)
        msg = f"Unknown tree generation mode {tree_mode!r}, expected one of {allowed}."
        raise ValueError(msg) from None


def random_tree(
    num_taxa: int,
    tree_mode: TreeGenMode | str,
    rand_seed: int | None = None,
) -> PhyloNode:
    """Generate a random phylogenetic tree.

    Generates a random tree through IQ-TREE.

    Parameters
    ----------
    num_taxa : int
        The number of taxa per tree.
    tree_mode : TreeGenMode | str
        How the tree is generated. Strings are the TreeGenMode member
        names (case-insensitive), e.g. "YULE_HARDING".
    rand_seed : int | None, optional
        The random seed - None means no seed is used, by default None.

    Returns
    -------
    PhyloNode
        A random phylogenetic tree.

    """
    tree_mode = _get_tree_gen_mode(tree_mode)
    rand_seed = process_rand_seed_nonzero(rand_seed)

    newick = iq_random_tree(
        num_taxa=num_taxa,
        tree_gen_mode=tree_mode.name,
        num_trees=1,
        rand_seed=rand_seed,
    ).strip()
    return make_tree(newick)
