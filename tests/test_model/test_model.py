import re

import pytest

from piqtree.model import (
    AaModel,
    CustomBaseFreq,
    DiscreteGammaModel,
    FreeRateModel,
    FreqType,
    LieModel,
    Model,
    RateModel,
    StandardDnaModel,
    SubstitutionModel,
    make_model,
)


def check_make_model(
    sub_mod: SubstitutionModel,
    freq_type: FreqType | CustomBaseFreq | None = None,
    rate_model: RateModel | None = None,
    *,
    invariable_sites: bool | float = False,
) -> None:
    model = Model(sub_mod, freq_type, rate_model, invariable_sites=invariable_sites)
    expected = str(model)

    # Check the expected string is approximately generated correctly

    # float and True should be True, False should be False
    assert (invariable_sites is not False) == model.invariable_sites

    if not isinstance(invariable_sites, bool):
        assert re.search(r"\+I\{\d+(\.\d+)?\}", expected)
        assert model.proportion_invariable_sites == invariable_sites
        assert model.invariable_sites
    else:
        assert model.proportion_invariable_sites is None
        assert model.invariable_sites == invariable_sites
        assert ("+I" in expected) == invariable_sites

    if isinstance(rate_model, DiscreteGammaModel):
        assert isinstance(model.rate_model, DiscreteGammaModel)
        assert "+G" in expected
        if rate_model.alpha is not None:
            assert re.search(r"\+G\d*\{\d+(\.\d+)?\}", expected)
    else:
        assert not isinstance(model.rate_model, DiscreteGammaModel)
        assert "+G" not in expected

    if isinstance(rate_model, FreeRateModel):
        assert isinstance(model.rate_model, FreeRateModel)
        assert "+R" in expected
        if rate_model.rates is not None:
            assert re.search(r"\+R\d*\{\d+(\.\d+)?(,\d+(\.\d+)?)*\}", expected)
    else:
        assert not isinstance(model.rate_model, FreeRateModel)
        assert "+R" not in expected

    if freq_type is not None:
        assert "+F" in expected
        if isinstance(freq_type, CustomBaseFreq):
            assert re.search(r"\+F\{\d+(\.\d+)?(,\d+(\.\d+)?)*\}", expected)
    else:
        assert "+F" not in expected

    # Check make_model
    got = str(make_model(expected))
    assert got == expected


@pytest.mark.parametrize(
    "sub_mod",
    [
        *StandardDnaModel.iter_available_models(),
        StandardDnaModel.GTR([4.39, 5.30, 4.39, 1.0, 12.1]),
        StandardDnaModel.TIM2([4.39, 5.30, 12.1]),
        *LieModel.iter_available_models(),
        *AaModel.iter_available_models(),
    ],
)
def test_make_model_sub_mod(sub_mod: SubstitutionModel) -> None:
    check_make_model(sub_mod)


@pytest.mark.parametrize(
    "sub_mod",
    [
        StandardDnaModel.iter_available_models()[0],
        LieModel.iter_available_models()[0],
        AaModel.iter_available_models()[0],
    ],
)
@pytest.mark.parametrize(
    "freq_type",
    [*list(FreqType), CustomBaseFreq([0.1, 0.2, 0.3, 0.4])],
)
def test_make_model_freq_type(
    sub_mod: SubstitutionModel,
    freq_type: FreqType | CustomBaseFreq,
) -> None:
    check_make_model(sub_mod, freq_type=freq_type)


@pytest.mark.parametrize(
    "sub_mod",
    [
        StandardDnaModel.iter_available_models()[-1],
        LieModel.iter_available_models()[-1],
        AaModel.iter_available_models()[-1],
    ],
)
@pytest.mark.parametrize("invariable_sites", [0.2, True])
def test_make_model_invar_sites(
    sub_mod: SubstitutionModel,
    invariable_sites: bool | float,
) -> None:
    check_make_model(sub_mod, invariable_sites=invariable_sites)


@pytest.mark.parametrize(
    "sub_mod",
    [
        StandardDnaModel.iter_available_models()[1],
        LieModel.iter_available_models()[1],
        AaModel.iter_available_models()[1],
    ],
)
@pytest.mark.parametrize(
    "rate_model",
    [
        DiscreteGammaModel(),
        DiscreteGammaModel(6),
        DiscreteGammaModel(2, 0.3),
        FreeRateModel(),
        FreeRateModel(6),
        FreeRateModel(2, [0.4, 0.6], [0.9, 0.1]),
    ],
)
def test_make_model_rate_model(
    sub_mod: SubstitutionModel,
    rate_model: RateModel,
) -> None:
    check_make_model(sub_mod, rate_model=rate_model)


@pytest.mark.parametrize(
    "sub_mod",
    [
        *StandardDnaModel.iter_available_models()[:3],
        StandardDnaModel.GTR([1.0, 2.0, 1.5, 3.7, 2.8]),
        *LieModel.iter_available_models()[:3],
        *AaModel.iter_available_models()[:3],
    ],
)
@pytest.mark.parametrize(
    "freq_type",
    [None, FreqType.FQ, CustomBaseFreq([0.1, 0.2, 0.3, 0.4])],
)
@pytest.mark.parametrize("invariable_sites", [False, 0.2, True])
@pytest.mark.parametrize(
    "rate_model",
    [
        None,
        DiscreteGammaModel(6),
        DiscreteGammaModel(2, 0.3),
        FreeRateModel(),
        FreeRateModel(2, [0.4, 0.6], [0.9, 0.1]),
    ],
)
def test_make_model(
    sub_mod: SubstitutionModel,
    freq_type: FreqType | CustomBaseFreq | None,
    invariable_sites: bool | float,
    rate_model: RateModel | None,
) -> None:
    check_make_model(sub_mod, freq_type, rate_model, invariable_sites=invariable_sites)


def test_model_repr() -> None:
    model = make_model("GTR")
    assert repr(model) == "Model(submod_type=GTR, freq_type=None, rate_type=None)"

    model = make_model("MK6.6+FO+R3")
    assert repr(model) == "Model(submod_type=MK6.6, freq_type=FO, rate_type=R3)"

    model = make_model("12.12+G")
    assert repr(model) == "Model(submod_type=12.12, freq_type=None, rate_type=G)"

    model = make_model("NQ.mammal+FQ")
    assert repr(model) == "Model(submod_type=NQ.mammal, freq_type=FQ, rate_type=None)"


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


def test_multiple_invariable_sites() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Model 'GTR+I+I' contains multiple specifications for invariable sites.",
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
