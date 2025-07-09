import re
from collections.abc import Sequence

import pytest

from piqtree.model import (
    AaModel,
    LieModel,
    Model,
    StandardDnaModel,
    StandardDnaModelInstance,
    SubstitutionModel,
    get_substitution_model,
)


@pytest.mark.parametrize("model_class", [StandardDnaModel, AaModel])
def test_number_of_descriptions(
    model_class: type[StandardDnaModel] | type[AaModel],
) -> None:
    assert len(model_class) == len(model_class._descriptions())


@pytest.mark.parametrize("model_class", [StandardDnaModel, AaModel])
def test_descriptions_exist(
    model_class: type[StandardDnaModel] | type[AaModel],
) -> None:
    for model in model_class:
        # Raises an error if description not present
        _ = model.description


@pytest.mark.parametrize(
    ("model_class", "model_type"),
    [(StandardDnaModel, "nucleotide"), (AaModel, "protein")],
)
def test_model_type(
    model_class: type[StandardDnaModel] | type[AaModel],
    model_type: str,
) -> None:
    assert model_class.model_type() == model_type

    for model in model_class:
        assert model.model_type() == model_type


@pytest.mark.parametrize(
    ("submod_type", "iqtree_str"),
    [
        (StandardDnaModel.F81, "F81"),
        (LieModel.LIE_10_34, "10.34"),
        (LieModel.LIE_4_4a("WS"), "WS4.4a"),
        (LieModel.LIE_5_11c("RY"), "RY5.11c"),
        (AaModel.NQ_insect, "NQ.insect"),
        ("NQ.yeast", "NQ.yeast"),
        ("GTR", "GTR"),
        ("2.2b", "2.2b"),
    ],
)
def test_get_substitution_model(
    submod_type: SubstitutionModel | str,
    iqtree_str: str,
) -> None:
    out = get_substitution_model(submod_type)
    assert isinstance(out, SubstitutionModel)
    assert out.iqtree_str() == iqtree_str


@pytest.mark.parametrize(
    ("model_str", "params"),
    [("GTR", (1.0, 2.0, 1.5, 3.7, 2.8)), ("TIM2", (4.39, 5.30, 12.1))],
)
def test_model_params(model_str: str, params: Sequence[float]) -> None:
    iqtree_str = model_str + f"{{{','.join(map(str, params))}}}"
    model = get_substitution_model(iqtree_str)

    assert model.iqtree_str() == iqtree_str
    assert isinstance(model, StandardDnaModelInstance)
    assert model.model_params is not None
    for model_param, param in zip(model.model_params, params, strict=True):
        assert model_param == param


def test_missing_bracket() -> None:
    model = "GTR{4.39,5.30,4.39,1.0,12.1"
    with pytest.raises(
        ValueError,
        match=re.escape(f"Missing closing bracket for parameterisation of {model!r}"),
    ):
        _ = get_substitution_model(model)


@pytest.mark.parametrize(
    "submod_type",
    ["FQ", "F", "+GTR", "AA", "G8", ""],
)
def test_invalid_substitution_model(submod_type: str) -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(f"Unknown substitution model: {submod_type!r}"),
    ):
        _ = get_substitution_model(submod_type)


def test_bad_lie_model() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("Invalid Lie Model pairing prefix: 'RS'"),
    ):
        _ = LieModel.LIE_10_34("RS")  # type: ignore[arg-type]


def test_lie_model_enum() -> None:
    # Test that using just the enum instead of LieModelInstance still works
    lie_model = LieModel.LIE_3_3a
    assert str(Model(lie_model)) == "3.3a"
    assert lie_model.model_type() == "nucleotide"
    assert LieModel.model_type() == "nucleotide"
    assert (
        lie_model.description
        == "Reversible model. Equal base frequencies. equiv. to K3P"
    )


@pytest.mark.parametrize(
    ("model", "base"),
    [
        (AaModel.WAG, AaModel.WAG),
        (LieModel.LIE_2_2b("RY"), LieModel.LIE_2_2b),
        (LieModel.LIE_10_12, LieModel.LIE_10_12),
        (StandardDnaModel.HKY([0.2]), StandardDnaModel.HKY),
        (StandardDnaModel.K3P, StandardDnaModel.K3P),
    ],
)
def test_base_model(model: SubstitutionModel, base: SubstitutionModel) -> None:
    assert model.base_model == base


def test_extended_descriptions() -> None:
    for model in StandardDnaModel:
        assert model.description == model().description


def test_extended_num_available() -> None:
    assert (
        StandardDnaModel.num_available_models()
        == StandardDnaModelInstance.num_available_models()
    )


def test_extended_iter_available() -> None:
    for m1, m2 in zip(
        StandardDnaModel.iter_available_models(),
        StandardDnaModelInstance.iter_available_models(),
        strict=True,
    ):
        assert m1 == m2


def test_extended_model_type() -> None:
    assert StandardDnaModelInstance.model_type() == "nucleotide"


def test_base_does_not_work() -> None:
    with pytest.raises(NotImplementedError):
        _ = SubstitutionModel.model_type()
    with pytest.raises(NotImplementedError):
        _ = SubstitutionModel().iqtree_str()
    with pytest.raises(NotImplementedError):
        _ = SubstitutionModel().base_model
    with pytest.raises(NotImplementedError):
        _ = SubstitutionModel.iter_available_models()
    with pytest.raises(NotImplementedError):
        _ = SubstitutionModel.num_available_models()
    with pytest.raises(NotImplementedError):
        _ = SubstitutionModel().description
