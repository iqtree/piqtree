import contextlib
import functools
from collections.abc import Sequence
from enum import Enum, unique


@unique
class FreqType(Enum):
    """Types of base frequencies."""

    F = "F"
    FO = "FO"
    FQ = "FQ"

    @staticmethod
    @functools.cache
    def _descriptions() -> dict["FreqType", str]:
        return {
            FreqType.F: "Empirical state frequency observed from the data.",
            FreqType.FO: "State frequency optimized by maximum-likelihood from the data. Note that this is with letter-O and not digit-0.",
            FreqType.FQ: "Equal state frequency.",
        }

    @property
    def description(self) -> str:
        """The description of the FreqType.

        Returns
        -------
        str
            The description of the FreqType.

        """
        return self._descriptions()[self]

    def iqtree_str(self) -> str:
        return self.value


class CustomBaseFreq:
    def __init__(self, frequencies: Sequence[float]) -> None:
        """Create a custom base frequency specification.

        For DNA models the order of the frequencies are:

        A C G T

        For AA models the order of the frequencies are:

        | A   | R   | N   | D   | C   | Q  | E  | G  | H  | I  | L   | K   | M   | F   | P   | S   | T   | W   | Y   | V   |
        |-----|-----|-----|-----|-----|----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
        | Ala | Arg | Asn | Asp | Cys | Gln| Glu| Gly| His| Ile| Leu | Lys | Met | Phe | Pro | Ser | Thr | Trp | Tyr | Val |

        IQ-TREE normalises the frequencies if they do not sum to one.

        Parameters
        ----------
        frequencies : Sequence[float]
            The fixed base frequencies.

        """
        self.frequencies = frequencies

        if len(self.frequencies) not in (4, 20):
            msg = f"Expected either 4 frequencies for DNA model or 20 for AA model but got {len(self.frequencies)}"
            raise ValueError(msg)

    @property
    def description(self) -> str:
        return "Fixed based frequencies."

    def iqtree_str(self) -> str:
        return f"F{{{','.join(map(str, self.frequencies))}}}"

    @classmethod
    def from_str(cls, base_freq_str: str) -> "CustomBaseFreq":
        if len(base_freq_str) == 0:
            msg = "An empty string does not specify base frequencies."
            raise ValueError(msg)

        if base_freq_str[0] != "F":
            msg = f"A CustomBaseFreq must start with F but got {base_freq_str!r}."
            raise ValueError(msg)

        if base_freq_str[1] != "{" or base_freq_str[-1] != "}":
            msg = f"A CustomBaseFreq must be parameterised with {{}} but got {base_freq_str!r}."
            raise ValueError(msg)

        try:
            frequencies = tuple(
                float(val.strip()) for val in base_freq_str[2:-1].split(",")
            )
        except ValueError:
            msg = f"Unable to parse parameters for CustomBaseFreq: {base_freq_str!r}."
            raise ValueError(msg) from None

        return cls(frequencies)


def get_freq_type(
    base_freq_str: str | FreqType | CustomBaseFreq,
) -> FreqType | CustomBaseFreq:
    """Return the FreqType enum for a given name.

    Parameters
    ----------
    base_freq_str : str | FreqType | CustomBaseFreq
        The string form of the frequency type.

    Returns
    -------
    FreqType | CustomBaseFreq
        The resolved FreqType Enum or CustomBaseFreq.

    Raises
    ------
    ValueError
        If the FreqType name cannot be resolved.

    """
    if isinstance(base_freq_str, (FreqType, CustomBaseFreq)):
        return base_freq_str

    shortened = base_freq_str.lstrip("+")

    with contextlib.suppress(KeyError):
        return FreqType[shortened]

    if shortened.startswith("F{"):
        return CustomBaseFreq.from_str(shortened)

    msg = f"Unknown state frequency type: {base_freq_str!r}."
    raise ValueError(msg)
