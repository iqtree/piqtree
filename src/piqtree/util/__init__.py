import secrets

from cogent3.core.tree import PhyloNode


def get_newick(tree: PhyloNode) -> str:
    return tree.get_newick(with_distances=True, escape_name=False)


def make_rand_seed() -> int:
    """Make a 32-bit random seed.

    Returns
    -------
    int
        A 32-bit random seed.
    """
    seed = secrets.randbits(32)
    return seed if seed < (1 << 31) else seed - (1 << 32)


def make_nonzero_rand_seed() -> int:
    """Make a non-zero 32-bit random seed.

    Returns
    -------
    int
        A non-zero 32-bit random seed.
    """
    while (seed := make_rand_seed()) == 0:
        pass  # pragma: no cover
    return seed


def process_rand_seed_nonzero(seed: int | None) -> int:
    """Process a random seed before sending to IQ-TREE.

    Some functions treat 0 as non-deterministic. For
    consistency between methods, 0 is treated as always
    deterministic and replaced with a pre-determined
    number.

    If the seed is None, a random seed is generated.
    If the seed is 0, it is replaced with a pre-determined
    number (1074213633).
    If the seed is anything else, returns that number.

    Parameters
    ----------
    seed : int | None
        The random seed to process.

    Returns
    -------
    int
        The original seed if no processing is required.
        Otherwise if the original seed was None, a random non-zero seed
        is returned. If the original seed was zero, a pre-determined
        number is returned.
    """
    if seed == 0:
        seed = 1074213633  # Randomly chosen once with make_rand_seed
    if seed is None:
        seed = make_nonzero_rand_seed()
    return seed
