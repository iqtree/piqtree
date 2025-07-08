import re

import pytest

from piqtree.model import (
    DiscreteGammaModel,
    FreeRateModel,
    RateModel,
    get_rate_type,
)


def test_rate_model_uninstantiable() -> None:
    with pytest.raises(TypeError):
        _ = RateModel()  # type: ignore[abstract]


@pytest.mark.parametrize(
    ("invariable_sites", "rate_model", "iqtree_str"),
    [
        (False, None, ""),
        (True, None, "I"),
        (0.1, None, "I{0.1}"),
        (False, DiscreteGammaModel(), "G"),
        (True, DiscreteGammaModel(), "I+G"),
        (False, FreeRateModel(), "R"),
        (True, FreeRateModel(), "I+R"),
        (False, DiscreteGammaModel(8), "G8"),
        (True, DiscreteGammaModel(8), "I+G8"),
        (False, FreeRateModel(8), "R8"),
        (True, FreeRateModel(8), "I+R8"),
        (False, DiscreteGammaModel(2, 0.5), "G2{0.5}"),
        (0.3, DiscreteGammaModel(3, 0.7), "I{0.3}+G3{0.7}"),
        (False, FreeRateModel(2, [0.3, 0.7], [0.1, 0.9]), "R2{0.3,0.1,0.7,0.9}"),
        (
            0.2,
            FreeRateModel(3, [0.2, 0.3, 0.5], [2, 0.1, 1.3]),
            "I{0.2}+R3{0.2,2,0.3,0.1,0.5,1.3}",
        ),
        (False, "G", "G"),
        (True, "+G", "I+G"),
        (False, "+R", "R"),
        (True, "R", "I+R"),
        (False, "G8", "G8"),
        (True, "+G8", "I+G8"),
        (False, "+R8", "R8"),
        (True, "R8", "I+R8"),
        (False, "+G42", "G42"),
        (True, "G42", "I+G42"),
        (False, "R42", "R42"),
        (True, "+R42", "I+R42"),
    ],
)
def test_get_rate_type(
    invariable_sites: bool | float,
    rate_model: RateModel | None,
    iqtree_str: str,
) -> None:
    model = get_rate_type(invariable_sites=invariable_sites, rate_model=rate_model)
    assert model.iqtree_str() == iqtree_str

    if rate_model is None:
        model = get_rate_type(invariable_sites=invariable_sites)
        assert model.iqtree_str() == iqtree_str

    if not invariable_sites:
        model = get_rate_type(rate_model=rate_model)
        assert model.iqtree_str() == iqtree_str


@pytest.mark.parametrize("invariable_sites", [True, False])
@pytest.mark.parametrize(
    "bad_rate_model",
    ["M", "T46"],
)
def test_invalid_rate_model_name(
    invariable_sites: bool,
    bad_rate_model: str,
) -> None:
    with pytest.raises(
        ValueError,
        match=f"Unexpected value for rate_model {bad_rate_model!r}",
    ):
        _ = get_rate_type(invariable_sites=invariable_sites, rate_model=bad_rate_model)


def test_bad_rate_categories() -> None:
    bad_rate_model = "R2D2"
    expected_error = re.escape(
        f"Invalid specification for rate categories {bad_rate_model!r}",
    )
    with pytest.raises(
        ValueError,
        match=expected_error,
    ):
        _ = get_rate_type(invariable_sites=True, rate_model=bad_rate_model)


@pytest.mark.parametrize("invariable_sites", [True, False])
@pytest.mark.parametrize(
    "bad_rate_model",
    [4, 3.15, ["R3", "G"]],
)
def test_invalid_rate_model_type(
    invariable_sites: bool,
    bad_rate_model: float | list,
) -> None:
    with pytest.raises(
        TypeError,
        match=f"Unexpected type for rate_model: {type(bad_rate_model)}",
    ):
        _ = get_rate_type(
            invariable_sites=invariable_sites,
            rate_model=bad_rate_model,  # type: ignore[arg-type]
        )


def test_discrete_empty_str() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("An empty string is not a DiscreteGammaModel."),
    ):
        _ = DiscreteGammaModel.from_str("")


def test_discrete_bad_start() -> None:
    rate_model_str = "R"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"A DiscreteGammaModel must start with G but got {rate_model_str!r}.",
        ),
    ):
        _ = DiscreteGammaModel.from_str(rate_model_str)


def test_discrete_missing_bracket() -> None:
    rate_model_str = "G{0.5"
    with pytest.raises(
        ValueError,
        match=re.escape(f"Missing end bracket for parameterisation {rate_model_str!r}"),
    ):
        _ = DiscreteGammaModel.from_str(rate_model_str)


def test_discrete_bad_parameter() -> None:
    rate_model_str = "G{cat}"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Parameterisation of Discrete Gamma Model is not a number {rate_model_str!r}",
        ),
    ):
        _ = DiscreteGammaModel.from_str(rate_model_str)


def test_free_empty_str() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("An empty string is not a FreeRateModel."),
    ):
        _ = FreeRateModel.from_str("")


def test_free_bad_start() -> None:
    rate_model_str = "G"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"A FreeRateModel must start with R but got {rate_model_str!r}.",
        ),
    ):
        _ = FreeRateModel.from_str(rate_model_str)


def test_free_missing_bracket() -> None:
    rate_model_str = "R2{0.5,1,0.5,1"
    with pytest.raises(
        ValueError,
        match=re.escape(f"Missing end bracket for parameterisation {rate_model_str!r}"),
    ):
        _ = FreeRateModel.from_str(rate_model_str)


def test_free_bad_parameter() -> None:
    rate_model_str = "R2{0.5,cat,0.5,1}"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Unable to parse parameters for FreeRateModel: {rate_model_str!r}.",
        ),
    ):
        _ = FreeRateModel.from_str(rate_model_str)


def test_free_bad_number_of_parameters() -> None:
    rate_model_str = "R2{0.5,1,0.5,1,0.5,1}"
    with pytest.raises(
        ValueError,
        match=re.escape("Expected 4 parameters but got 6."),
    ):
        _ = FreeRateModel.from_str(rate_model_str)

    rate_model_str = "R2{0.5,1,0.5}"
    with pytest.raises(
        ValueError,
        match=re.escape("Expected 4 parameters but got 3."),
    ):
        _ = FreeRateModel.from_str(rate_model_str)
