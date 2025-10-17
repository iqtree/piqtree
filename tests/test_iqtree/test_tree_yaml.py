import re
from typing import Any

import pytest
from cogent3 import make_tree

from piqtree import make_model
from piqtree.exceptions import ParseIqTreeError
from piqtree.iqtree._tree import _process_tree_yaml, _tree_equal


@pytest.fixture
def newick_not_in_candidates() -> list[dict[str, Any]]:
    return [
        {  # Newick string not in candidate set (different branch length)
            "CandidateSet": {
                0: "-6519.33018689 (0:0.0058955371,1:0.0026486308,(2:0.0230933557,3:0.3069062230):0.01387802789);",
                1: "-6540.1924365 (0:0.0059276645,(1:0.0026255655,2:0.0369876991):2.43436209e-06,3:0.3205542282);",
                2: "-6540.32542968 (0:0.0059468612,(1:0.0021841363,3:0.3203544844):2.076530752e-06,2:0.0369270512);",
            },
            "ModelDNA": {
                "rates": "1, 3.815110072, 1, 1, 3.815110072, 1",
                "state_freq": "0.3640205807, 0.1862366777, 0.217291437, 0.2324513047",
            },
            "PhyloTree": {
                "newick": "(0:0.0068955371,1:0.0026486308,(2:0.0230933557,3:0.3069062230):0.01387802789);",
            },
            "StopRule": {
                "curIteration": 101,
                "start_real_time": 1731027583,
                "time_vec": None,
            },
            "boot_consense_logl": 0,
            "contree_rfdist": -1,
            "finished": True,
            "finishedCandidateSet": True,
            "finishedModelFinal": True,
            "finishedModelInit": True,
            "initTree": "(0:0.0063001853,1:0.0022115739,(2:0.0203850510,3:0.3497395366):0.01883712169);",
            "iqtree": {
                "seed": 598834595,
                "start_time": 1731027582,
                "version": "2.3.6.lib",
            },
        },
        {  # Newick string not in candidate set (different names)
            "CandidateSet": {
                0: "-6519.33018689 (0:0.0058955371,2:0.0026486308,(1:0.0230933557,3:0.3069062230):0.01387802789);",
                1: "-6540.1924365 (0:0.0059276645,(1:0.0026255655,2:0.0369876991):2.43436209e-06,3:0.3205542282);",
                2: "-6540.32542968 (0:0.0059468612,(1:0.0021841363,3:0.3203544844):2.076530752e-06,2:0.0369270512);",
            },
            "ModelDNA": {
                "rates": "1, 3.815110072, 1, 1, 3.815110072, 1",
                "state_freq": "0.3640205807, 0.1862366777, 0.217291437, 0.2324513047",
            },
            "PhyloTree": {
                "newick": "(0:0.0058955371,1:0.0026486308,(2:0.0230933557,3:0.3069062230):0.01387802789);",
            },
            "StopRule": {
                "curIteration": 101,
                "start_real_time": 1731027583,
                "time_vec": None,
            },
            "boot_consense_logl": 0,
            "contree_rfdist": -1,
            "finished": True,
            "finishedCandidateSet": True,
            "finishedModelFinal": True,
            "finishedModelInit": True,
            "initTree": "(0:0.0063001853,1:0.0022115739,(2:0.0203850510,3:0.3497395366):0.01883712169);",
            "iqtree": {
                "seed": 598834595,
                "start_time": 1731027582,
                "version": "2.3.6.lib",
            },
        },
    ]


@pytest.fixture
def non_lie_dna_with_rate_model() -> dict[str, Any]:
    return {
        "CandidateSet": {
            0: "-6736.94578464 (0:0.0063211201,1:0.0029675780,(2:0.0228519739,3:0.3072009029):0.01373649616);",
            1: "-6757.78815651 (0:0.0063607954,(1:0.0030079874,2:0.0365597715):2.296825575e-06,3:0.3208135518);",
            2: "-6758.07765021 (0:0.0063826033,(1:0.0021953253,3:0.3207201830):0.0001145372551,2:0.0365362763);",
        },
        "ModelDNA": {
            "rates": "1, 3.82025079, 1, 1, 3.82025079, 1",
            "state_freq": "0.3628523161, 0.1852938562, 0.2173913044, 0.2344625233",
        },
        "PhyloTree": {
            "newick": "(0:0.0063211201,1:0.0029675780,(2:0.0228519739,3:0.3072009029):0.01373649616);",
        },
        "RateGammaInvar": {"gamma_shape": 1.698497993, "p_invar": 1.002841144e-06},
        "StopRule": {
            "curIteration": 101,
            "start_real_time": 1724397157,
            "time_vec": None,
        },
        "boot_consense_logl": 0,
        "contree_rfdist": -1,
        "finished": True,
        "finishedCandidateSet": True,
        "finishedModelFinal": True,
        "finishedModelInit": True,
        "initTree": "(0:0.0063680036,1:0.0026681490,(2:0.0183861083,3:0.3034074996):0.01838610827);",
        "iqtree": {"seed": 95633264, "start_time": 1724397157, "version": "2.3.6.lib"},
    }


@pytest.fixture
def lie_dna_model() -> dict[str, Any]:
    return {
        "CandidateSet": {
            0: "-6606.82337571 ((0:0.0058632731,(2:0.0225450645,3:0.3055011501):0.01414755595):2.340509431e-06,1:0.0026433577):0;",
            1: "-6606.82337581 (0:0.0058609304,(1:0.0026456850,(2:0.0225450645,3:0.3055011501):0.01414755618):2.340509431e-06):0;",
            2: "-6606.82337581 ((0:0.0058632731,1:0.0026456850):2.340509431e-06,(2:0.0225450645,3:0.3055011501):0.01414520892):0;",
            3: "-6606.82337594 (((0:0.0058630055,1:0.0026461740):0.01414791335,2:0.0225451031):0.1178978025,3:0.1876035160):0;",
            4: "-6628.65519886 (0:0.0000023102,((1:0.0027556399,3:0.3193110364):2.07651076e-06,2:0.0366286909):0.005770021808):1e-06;",
            5: "-6628.70596445 (0:0.0000022647,((1:0.0027641463,2:0.0366284785):2.43434433e-06,3:0.3193357780):0.005401508376):1e-06;",
        },
        "ModelLieMarkovRY2.2b": {
            "model_parameters": 0.4841804549,
            "state_freq": "0.25, 0.25, 0.25, 0.25",
        },
        "PhyloTree": {
            "newick": "((0:0.0058632731,(2:0.0225450645,3:0.3055011501):0.01414755595):2.340509431e-06,1:0.0026433577):0;",
        },
        "StopRule": {
            "curIteration": 101,
            "start_real_time": 1729226726,
            "time_vec": None,
        },
        "boot_consense_logl": 0,
        "contree_rfdist": -1,
        "finished": True,
        "finishedCandidateSet": True,
        "finishedModelFinal": True,
        "finishedModelInit": True,
        "initTree": "(0:0.0063001853,1:0.0022115739,(2:0.0203850510,3:0.3497395366):0.01883712169);",
        "iqtree": {"seed": 1, "start_time": 1729226725, "version": "2.3.6.lib"},
    }


@pytest.fixture
def unrest_model() -> dict[str, Any]:
    return {
        "CandidateSet": {
            0: "-6736.94578464 (0:0.0063211201,1:0.0029675780,(2:0.0228519739,3:0.3072009029):0.01373649616);",
            1: "-6757.78815651 (0:0.0063607954,(1:0.0030079874,2:0.0365597715):2.296825575e-06,3:0.3208135518);",
            2: "-6758.07765021 (0:0.0063826033,(1:0.0021953253,3:0.3207201830):0.0001145372551,2:0.0365362763);",
        },
        "ModelUnrest": {
            "rates": "0.9543854572, 2.942336813, 0.5526681694, 0.7599462192, 0.897052927, 3.38882313, 4.630689862, 0.7223043133, 0.4970030179, 0.8401385788, 2.854877074, 1",
            "state_freq": "0.3234606464, 0.2202693695, 0.2343932981, 0.221876686",
        },
        "PhyloTree": {
            "newick": "(0:0.0063211201,1:0.0029675780,(2:0.0228519739,3:0.3072009029):0.01373649616);",
        },
    }


def test_newick_not_in_candidates(
    newick_not_in_candidates: list[dict[str, Any]],
) -> None:
    for yaml in newick_not_in_candidates:
        with pytest.raises(
            ParseIqTreeError,
            match=re.escape("IQ-TREE output malformated, likelihood not found."),
        ):
            _ = _process_tree_yaml(yaml, ["a", "b", "c", "d"], make_model("JC"))


def test_non_lie_dna_with_rate_model(
    non_lie_dna_with_rate_model: dict[str, Any],
) -> None:
    # test parsing yaml file containing fits from non-Lie DnaModel and rate heterogeneity model
    edge_params = {
        "rates": {
            "A/C": 1.0,
            "A/G": 3.82025079,
            "A/T": 1.0,
            "C/G": 1.0,
            "C/T": 3.82025079,
        },
        "mprobs": {
            "A": 0.3628523161,
            "C": 0.1852938562,
            "G": 0.2173913044,
            "T": 0.2344625233,
        },
    }
    rate_params = {"gamma_shape": 1.698497993, "p_invar": 1.002841144e-06}
    tree = _process_tree_yaml(
        non_lie_dna_with_rate_model,
        ["a", "b", "c", "d"],
        make_model("GTR+I+G"),
    )
    print(tree[0].params)
    for rate, value in edge_params["rates"].items():
        assert tree[0].params[rate] == value
    assert tree.params["RateGammaInvar"] == rate_params


def test_non_lie_dna_model_motif_absent(
    non_lie_dna_with_rate_model: dict[str, Any],
) -> None:
    non_lie_dna_with_rate_model["ModelDNA"].pop("state_freq")
    with pytest.raises(
        ParseIqTreeError,
        match=re.escape("IQ-TREE output malformated, motif parameters not found."),
    ):
        _ = _process_tree_yaml(
            non_lie_dna_with_rate_model,
            ["a", "b", "c", "d"],
            make_model("GTR+I+G"),
        )


def test_non_lie_dna_model_rate_absent(
    non_lie_dna_with_rate_model: dict[str, Any],
) -> None:
    non_lie_dna_with_rate_model["ModelDNA"].pop("rates")
    with pytest.raises(
        ParseIqTreeError,
        match=re.escape("IQ-TREE output malformated, rate parameters not found."),
    ):
        _ = _process_tree_yaml(
            non_lie_dna_with_rate_model,
            ["a", "b", "c", "d"],
            make_model("GTR+I+G"),
        )


def test_lie_dna_model(
    lie_dna_model: dict[str, Any],
) -> None:
    # test parsing yaml file containing fits from Lie DnaModel
    model_parameters = {
        "model_parameters": 0.4841804549,
        "mprobs": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    }
    tree = _process_tree_yaml(lie_dna_model, ["a", "b", "c", "d"], make_model("RY2.2b"))
    assert tree.params["ModelLieMarkovRY2.2b"] == model_parameters


def test_lie_dna_model_motif_absent(
    lie_dna_model: dict[str, Any],
) -> None:
    lie_dna_model["ModelLieMarkovRY2.2b"].pop("state_freq")
    with pytest.raises(
        ParseIqTreeError,
        match=re.escape("IQ-TREE output malformated, motif parameters not found."),
    ):
        _ = _process_tree_yaml(
            lie_dna_model,
            ["a", "b", "c", "d"],
            make_model("RY2.2b"),
        )


def test_unrest_model(
    unrest_model: dict[str, Any],
) -> None:
    tree = _process_tree_yaml(unrest_model, ["a", "b", "c", "d"], make_model("UNREST"))
    assert "A/C" in tree[0].params
    assert tree.params["mprobs"]


@pytest.mark.parametrize(
    ("candidate", "expected"),
    [
        ("((a:1.0,b:0.9),c:0.8);", True),
        ("((a:0.9,b:0.9),c:0.8);", False),
        ("((a:1.0,c:0.8),b:0.9);", False),
    ],
)
def test_tree_equal(candidate: str, expected: bool) -> None:
    tree = make_tree("((a:1.0,b:0.9),c:0.8);")
    candidate = make_tree(candidate)
    assert _tree_equal(tree, candidate) == expected
