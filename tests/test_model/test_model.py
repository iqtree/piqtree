import re

import pytest

from piqtree.model import (
    AaModel,
    DiscreteGammaModel,
    DnaModel,
    FreeRateModel,
    FreqType,
    Model,
    RateModel,
    SubstitutionModel,
    make_model,
)


@pytest.mark.parametrize("sub_mod", list(DnaModel) + list(AaModel))
@pytest.mark.parametrize("freq_type", [None, *list(FreqType)])
@pytest.mark.parametrize("invariant_sites", [False, True])
@pytest.mark.parametrize(
    "rate_model",
    [
        None,
        DiscreteGammaModel(),
        FreeRateModel(),
        DiscreteGammaModel(6),
        FreeRateModel(6),
    ],
)
def test_make_model(
    sub_mod: SubstitutionModel,
    freq_type: FreqType,
    invariant_sites: bool,
    rate_model: RateModel,
) -> None:
    model = Model(sub_mod, freq_type, rate_model, invariant_sites=invariant_sites)
    expected = str(model)

    # Check the expected string is approximately generated correctly
    if invariant_sites:
        assert "+I" in expected
    else:
        assert "+I" not in expected

    if isinstance(rate_model, DiscreteGammaModel):
        assert "+G" in expected
    else:
        assert "+G" not in expected

    if isinstance(rate_model, FreeRateModel):
        assert "+R" in expected
    else:
        assert "+R" not in expected

    if freq_type is not None:
        assert "+F" in expected
    else:
        assert "+F" not in expected

    # Check make_model
    got = str(make_model(expected))
    assert got == expected


def test_bad_sub_model() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("Unknown substitution model: 'GYR'"),
    ):
        make_model("GYR")


def test_multiple_freq_type() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Model 'GTR+FO+FO' contains multiple base frequency specifications.",
        ),
    ):
        make_model("GTR+FO+FO")


def test_multiple_invariant_sites() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Model 'GTR+I+I' contains multiple specifications for invariant sites.",
        ),
    ):
        make_model("GTR+I+I")


def test_multiple_rate_het() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Model 'GTR+G+R' contains multiple rate heterogeneity specifications.",
        ),
    ):
        make_model("GTR+G+R")


def test_unexpected_component() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Model 'GTR+Z' contains unexpected component.",
        ),
    ):
        make_model("GTR+Z")
