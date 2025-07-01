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
        (False, DiscreteGammaModel(), "G"),
        (True, DiscreteGammaModel(), "I+G"),
        (False, FreeRateModel(), "R"),
        (True, FreeRateModel(), "I+R"),
        (False, DiscreteGammaModel(8), "G8"),
        (True, DiscreteGammaModel(8), "I+G8"),
        (False, FreeRateModel(8), "R8"),
        (True, FreeRateModel(8), "I+R8"),
        (False, DiscreteGammaModel(42), "G42"),
        (True, DiscreteGammaModel(42), "I+G42"),
        (False, FreeRateModel(42), "R42"),
        (True, FreeRateModel(42), "I+R42"),
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
    invariable_sites: bool,
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
    ["M", "T46", "R2D2"],
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
