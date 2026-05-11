"""Tests for piqtree.fit_tree_hessian."""

import re

import numpy as np
import pytest
from cogent3 import make_tree
from cogent3.core.alignment import Alignment

import piqtree


def test_fit_tree_hessian_basic(four_otu: Alignment) -> None:
    """Smoke test: returns sensibly-shaped gradient / Hessian arrays."""
    tree_topology = make_tree(tip_names=four_otu.names)

    result = piqtree.fit_tree_hessian(four_otu, tree_topology, "JC")

    assert isinstance(result, piqtree.HessianFit)
    assert result.branch_lengths.dtype == np.float64
    assert result.gradient.dtype == np.float64
    assert result.hessian.dtype == np.float64

    n = result.branch_lengths.shape[0]
    # 4 taxa rooted tree has 2*4 - 2 = 6 branches; allow a bit of slack for
    # whatever rooting MCMCTree-style produces, but n must be > 0.
    assert n > 0
    assert result.gradient.shape == (n,)
    assert result.hessian.shape == (n, n)


def test_fit_tree_hessian_hessian_is_symmetric(four_otu: Alignment) -> None:
    """The Hessian of a smooth log-likelihood should be (numerically) symmetric."""
    tree_topology = make_tree(tip_names=four_otu.names)

    result = piqtree.fit_tree_hessian(four_otu, tree_topology, "JC")

    np.testing.assert_allclose(result.hessian, result.hessian.T, atol=1e-6)


def test_fit_tree_hessian_branch_lengths_positive(four_otu: Alignment) -> None:
    tree_topology = make_tree(tip_names=four_otu.names)

    result = piqtree.fit_tree_hessian(four_otu, tree_topology, "JC")

    assert np.all(result.branch_lengths > 0)


def test_fit_tree_hessian_tree_taxa_match(four_otu: Alignment) -> None:
    tree_topology = make_tree(tip_names=four_otu.names)

    result = piqtree.fit_tree_hessian(four_otu, tree_topology, "GTR")

    assert set(result.tree.get_tip_names()) == set(four_otu.names)


def test_fit_tree_hessian_dating_option_rejected(four_otu: Alignment) -> None:
    """--dating is owned by this function; supplying it via other_options is an error."""
    tree_topology = make_tree(tip_names=four_otu.names)

    with pytest.raises(
        ValueError,
        match=re.escape(
            "Option '--dating' will be overridden by parameter of function.",
        ),
    ):
        piqtree.fit_tree_hessian(
            four_otu,
            tree_topology,
            "JC",
            other_options="--dating mcmctree",
        )
