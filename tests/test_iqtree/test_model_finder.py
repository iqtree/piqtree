import multiprocessing

import pytest
from cogent3.core.alignment import Alignment

from piqtree.iqtree import ModelFinderResult, ModelResultValue, model_finder


def test_model_result_value_from_string() -> None:
    val = "123.45 10 0.678"
    result = ModelResultValue.from_string(val)
    assert result.lnL == 123.45
    assert result.nfp == 10
    assert result.tree_length == 0.678


@pytest.mark.parametrize(
    "bad_val",
    [
        "123.45.00 10 0.678",  # invalid float
        "123.45 10 10 0.678",  # too many values
        "123.45 10",  # too few values
    ],
)
def test_model_result_value_from_bad_string(bad_val: str) -> None:
    with pytest.raises(ValueError, match="Error parsing string"):
        _ = ModelResultValue.from_string(bad_val)


@pytest.mark.parametrize("model", ["GTR+F", "GTR+F+G"])
def test_model_finder_result(model: str) -> None:
    raw_data = {
        model: "123.45 10 0.678",
        "best_model_AIC": model,
        "best_model_AICc": model,
        "best_model_BIC": model,
        "best_tree_AIC": "((a,b),(c,d));",  # ignored
        "initTree": "((a,b),(c,d));",  # ignored
        "partition_type": 0,  # ignored
    }

    result = ModelFinderResult("test", raw_data)

    assert isinstance(result.model_stats[model], ModelResultValue)
    assert result.model_stats[model].lnL == 123.45
    assert result.model_stats[model].nfp == 10
    assert result.model_stats[model].tree_length == 0.678


def test_model_finder(five_otu: Alignment) -> None:
    result1 = model_finder(five_otu, rand_seed=1)
    result2 = model_finder(
        five_otu,
        num_threads=multiprocessing.cpu_count(),
        rand_seed=1,
    )
    assert str(result1.best_aic) == str(result2.best_aic)
    assert str(result1.best_aicc) == str(result2.best_aicc)
    assert str(result1.best_bic) == str(result2.best_bic)


def test_model_finder_restricted_submod(five_otu: Alignment) -> None:
    result = model_finder(five_otu, rand_seed=1, model_set={"HKY", "TIM"})
    assert str(result.best_aic).startswith(("HKY", "TIM"))
    assert str(result.best_aicc).startswith(("HKY", "TIM"))
    assert str(result.best_bic).startswith(("HKY", "TIM"))
