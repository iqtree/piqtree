from cogent3.core.tree import PhyloNode


def get_newick(tree: PhyloNode) -> str:
    return tree.get_newick(with_distances=True, escape_name=False)
