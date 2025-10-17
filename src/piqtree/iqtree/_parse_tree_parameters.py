from typing import Any

from cogent3.core.tree import PhyloNode

from piqtree.exceptions import ParseIqTreeError
from piqtree.model import Model, StandardDnaModel

# the order defined in IQ-TREE
# assume UNREST model has 12 rates, GTR and simpler models always have 6 rates present
RATE_PARS = "A/C", "A/G", "A/T", "C/G", "C/T", "G/T"
RATE_PARS_UNREST = (
    "A/C",
    "A/G",
    "A/T",
    "C/A",
    "C/G",
    "C/T",
    "G/A",
    "G/C",
    "G/T",
    "T/A",
    "T/C",
    "T/G",
)
MOTIF_PARS = "A", "C", "G", "T"


def _insert_edge_pars(tree: PhyloNode, **kwargs: dict) -> None:
    # inserts the edge parameters into each edge to match the structure of
    # PhyloNode
    for node in tree.get_edge_vector():
        # skip the rate parameters when node is the root
        if node.is_root():
            kwargs = {k: v for k, v in kwargs.items() if k == "mprobs"}
            del node.params["edge_pars"]
        node.params.update(kwargs)


def _edge_pars_for_cogent3(tree: PhyloNode, model: Model) -> None:
    base_model = model.submod_type.base_model

    rate_pars = tree.params["edge_pars"]["rates"]
    motif_pars = {"mprobs": tree.params["edge_pars"]["mprobs"]}
    # renames parameters to conform to cogent3's naming conventions
    if base_model in {StandardDnaModel.JC, StandardDnaModel.F81}:
        # skip rate_pars since rate parameters are constant in JC and F81
        _insert_edge_pars(
            tree,
            **motif_pars,
        )
        return
    if base_model in {StandardDnaModel.K80, StandardDnaModel.HKY}:
        rate_pars = {"kappa": rate_pars["A/G"]}

    elif base_model is StandardDnaModel.TN:
        rate_pars = {"kappa_r": rate_pars["A/G"], "kappa_y": rate_pars["C/T"]}

    elif base_model is StandardDnaModel.GTR:
        del rate_pars["G/T"]

    # applies global rate parameters to each edge
    _insert_edge_pars(
        tree,
        **rate_pars,
        **motif_pars,
    )


def _parse_nonlie_model(tree: PhyloNode, tree_yaml: dict) -> None:
    # parse motif and rate parameters in the tree_yaml for non-Lie DnaModel
    model_fits = tree_yaml.get("ModelDNA", {})

    state_freq_str = model_fits.get("state_freq", "")
    rate_str = model_fits.get("rates", "")

    # parse motif parameters, assign each to a name, and raise an error if not found
    if state_freq_str:
        # converts the strings of motif parameters into dictionary
        state_freq_list = [
            float(value) for value in state_freq_str.replace(" ", "").split(",")
        ]
        tree.params["edge_pars"] = {
            "mprobs": dict(zip(MOTIF_PARS, state_freq_list, strict=True)),
        }
    else:
        msg = "IQ-TREE output malformated, motif parameters not found."
        raise ParseIqTreeError(msg)

    # parse rate parameters, assign each to a name, and raise an error if not found
    if rate_str:
        rate_list = [float(value) for value in rate_str.replace(" ", "").split(",")]
        tree.params["edge_pars"]["rates"] = dict(
            zip(RATE_PARS, rate_list, strict=True),
        )
    else:
        msg = "IQ-TREE output malformated, rate parameters not found."
        raise ParseIqTreeError(msg)


def _parse_lie_model(
    tree: PhyloNode,
    tree_yaml: dict,
    lie_model_name: str,
) -> None:
    # parse motif and rate parameters in the tree_yaml for Lie DnaModel
    model_fits = tree_yaml.get(lie_model_name, {})

    # parse motif parameters, assign each to a name, and raise an error if not found
    state_freq_str = model_fits.get("state_freq", "")
    if state_freq_str:
        state_freq_list = [
            float(value) for value in state_freq_str.replace(" ", "").split(",")
        ]
        tree.params[lie_model_name] = {
            "mprobs": dict(zip(MOTIF_PARS, state_freq_list, strict=True)),
        }
    else:
        msg = "IQ-TREE output malformated, motif parameters not found."
        raise ParseIqTreeError(msg)

    # parse rate parameters, skipping LIE_1_1 (aka JC69) since its rate parameter is constant thus absent
    if "model_parameters" in model_fits:
        model_parameters = model_fits["model_parameters"]

        # convert model parameters to a list of floats if they are a string
        if isinstance(model_parameters, str):
            tree.params[lie_model_name]["model_parameters"] = [
                float(value) for value in model_parameters.replace(" ", "").split(",")
            ]
        else:
            # directly use the float
            tree.params[lie_model_name]["model_parameters"] = model_parameters


def _parse_unrest_model(tree: PhyloNode, tree_yaml: dict) -> None:
    model_fits = tree_yaml.get("ModelUnrest", {})

    state_freq_str = model_fits.get("state_freq", "")
    rate_str = model_fits.get("rates", "")

    # parse state frequencies
    if state_freq_str:
        state_freq_list = [
            float(value) for value in state_freq_str.replace(" ", "").split(",")
        ]
        tree.params["edge_pars"] = {
            "mprobs": dict(zip(MOTIF_PARS, state_freq_list, strict=True)),
        }
    else:
        msg = "IQ-TREE output malformated, motif parameters not found."
        raise ParseIqTreeError(msg)

    # parse rates
    if rate_str:
        rate_list = [float(value) for value in rate_str.replace(" ", "").split(",")]
        tree.params["edge_pars"]["rates"] = dict(
            zip(RATE_PARS_UNREST, rate_list, strict=True),
        )
    else:
        msg = "IQ-TREE output malformated, rate parameters not found."
        raise ParseIqTreeError(msg)


def parse_model_parameters(
    tree: PhyloNode,
    tree_yaml: dict[str, Any],
    model: Model,
) -> None:
    """Parse model parameters from the returned yaml format.

    Parameters
    ----------
    tree : PhyloNode
        The tree to attach model parameters to.
    tree_yaml : dict[str, Any]
        The yaml result returned from IQ-TREE.
    """
    # parse non-Lie DnaModel parameters
    if "ModelDNA" in tree_yaml:
        _parse_nonlie_model(tree, tree_yaml)

    elif "ModelUnrest" in tree_yaml:
        _parse_unrest_model(tree, tree_yaml)

    # parse Lie DnaModel parameters, handling various Lie model names
    elif key := next(
        (key for key in tree_yaml if key.startswith("ModelLieMarkov")),
        None,
    ):
        _parse_lie_model(tree, tree_yaml, key)

    # for non-Lie models, populate parameters to each branch and
    # rename them to mimic PhyloNode
    if "edge_pars" in tree.params:
        _edge_pars_for_cogent3(tree, model)
