import re
from collections.abc import Sequence

import pytest

from piqtree.model import CustomBaseFreq, FreqType, get_freq_type


def test_number_of_descriptions() -> None:
    assert len(FreqType) == len(FreqType._descriptions())


def test_descriptions_exist() -> None:
    for freq_type in FreqType:
        # Raises an error if description not present
        _ = freq_type.description

    _ = get_freq_type("F{0.1,0.2,0.3,0.4}").description


@pytest.mark.parametrize(
    ("freq_type", "iqtree_str"),
    [
        (FreqType.F, "F"),
        (FreqType.FO, "FO"),
        (FreqType.FQ, "FQ"),
        ("F", "F"),
        ("FO", "FO"),
        ("FQ", "FQ"),
        ("+F", "F"),
        ("+FO", "FO"),
        ("+FQ", "FQ"),
    ],
)
def test_get_freq_type(freq_type: FreqType | str, iqtree_str: str) -> None:
    out = get_freq_type(freq_type)
    assert isinstance(out, FreqType)
    assert out.iqtree_str() == iqtree_str


@pytest.mark.parametrize(
    ("freq_type", "iqtree_str"),
    [
        ("+F{0.1,0.2,0.3,0.4}", "F{0.1,0.2,0.3,0.4}"),
        (
            f"+F{{{','.join(map(str, [0.05] * 20))}}}",
            f"F{{{','.join(map(str, [0.05] * 20))}}}",
        ),
    ],
)
def test_get_custom_freq_type(freq_type: str, iqtree_str: str) -> None:
    out = get_freq_type(freq_type)
    assert isinstance(out, CustomBaseFreq)
    assert out.iqtree_str() == iqtree_str


@pytest.mark.parametrize(
    "freq_type",
    ["F0", "+F0", "+G", "+R9"],
)
def test_invalid_freq_type_name(freq_type: str) -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(f"Unknown state frequency type: {freq_type!r}"),
    ):
        _ = get_freq_type(freq_type)


@pytest.mark.parametrize(
    "freqs",
    [[0.1, 0.2, 0.7], [0.2] * 5, [0.04] * 25],
)
def test_invalid_freq_type_params(freqs: Sequence[float]) -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Expected either 4 frequencies for DNA model or 20 for AA model but got {len(freqs)}",
        ),
    ):
        _ = get_freq_type(f"F{{{','.join(map(str, freqs))}}}")


def test_custom_empty_str() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("An empty string does not specify base frequencies."),
    ):
        _ = CustomBaseFreq.from_str("")


def test_custom_bad_start() -> None:
    base_freq_str = "G{0.1,0.2,0.3,0.4}"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"A CustomBaseFreq must start with F but got {base_freq_str!r}.",
        ),
    ):
        _ = CustomBaseFreq.from_str(base_freq_str)


def test_custom_missing_bracket() -> None:
    base_freq_str = "F{0.1,0.2,0.3,0.4"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"A CustomBaseFreq must be parameterised with {{}} but got {base_freq_str!r}.",
        ),
    ):
        _ = CustomBaseFreq.from_str(base_freq_str)


def test_custom_bad_parameter() -> None:
    base_freq_str = "F{0.1,0.2,0.3,cat}"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Unable to parse parameters for CustomBaseFreq: {base_freq_str!r}.",
        ),
    ):
        _ = CustomBaseFreq.from_str(base_freq_str)
