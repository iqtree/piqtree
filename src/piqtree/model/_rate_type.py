from abc import ABC, abstractmethod


class RateModel(ABC):
    """Base class for rate models."""

    @abstractmethod
    def iqtree_str(self) -> str:
        """Convert to an iqtree settings string.

        Returns
        -------
        str
            String parsable by IQ-TREE for the rate heterogeneity model.

        """


class RateType:
    def __init__(
        self,
        *,
        invariable_sites: bool | float = False,
        rate_model: RateModel | None = None,
    ) -> None:
        """Rate heterogeneity across sites model.

        Parameters
        ----------
        invariable_sites : bool | float, optional
            Invariable sites, by default False.
            If a float in range [0,1) specifies the proportion of invariable sites.
        rate_model : RateModel | None, optional
            Discrete Gamma Model or FreeRate Model, by default None.

        """

        if not isinstance(invariable_sites, bool):
            if not (0 <= invariable_sites < 1):
                msg = "The proportion of invaraint sites must be in the range [0,1)"
                raise ValueError(msg)
            self.invariable_sites = True
            self.proportion_invariable: float | None = invariable_sites
        else:
            self.invariable_sites = invariable_sites
            self.proportion_invariable = None

        self.rate_model = rate_model

    def iqtree_str(self) -> str:
        """Convert to an iqtree settings string.

        Returns
        -------
        str
            String parsable by IQ-TREE for the rate heterogeneity model.

        """
        rate_type_str = "I" if self.invariable_sites else ""
        if self.proportion_invariable is not None:
            rate_type_str += f"{{{self.proportion_invariable}}}"

        if self.rate_model is None:
            return rate_type_str
        # Invariant sites and model need to be joined by a '+'
        if self.invariable_sites:
            rate_type_str += "+"
        return rate_type_str + self.rate_model.iqtree_str()

    @property
    def name(self) -> str:
        return self.iqtree_str()


class DiscreteGammaModel(RateModel):
    def __init__(
        self,
        rate_categories: int | None = None,
        alpha: float | None = None,
    ) -> None:
        """Discrete Gamma Model.

        Parameters
        ----------
        rate_categories : int, optional
            The number of rate categories, by default 4.
        alpha : float | None, optional
            The Gamma shape parameter - fixed if not None, by default None

        References
        ----------
        .. [1] Yang, Ziheng. "Maximum likelihood phylogenetic estimation from
           DNA sequences with variable rates over sites: approximate methods."
           Journal of Molecular evolution 39 (1994): 306-314.

        """
        self.rate_categories = rate_categories
        self.alpha = alpha

    def iqtree_str(self) -> str:
        out_str = "G"

        if self.rate_categories is not None:
            out_str += str(self.rate_categories)

        if self.alpha is not None:
            out_str += f"{{{self.alpha}}}"
        return out_str

    @classmethod
    def from_str(cls, rate_model_str: str) -> "DiscreteGammaModel":
        if len(rate_model_str) == 0:
            msg = "An empty string is not a DiscreteGammaModel."
            raise ValueError(msg)

        if rate_model_str[0] != "G":
            msg = f"A DiscreteGammaModel must start with G but got {rate_model_str!r}."
            raise ValueError(msg)

        if len(rate_model_str) == 1:
            return cls()

        rate_categories = _parse_rate_categories(rate_model_str)
        parameter_start = rate_model_str.find("{")

        # If there is no parameterisation
        if parameter_start == -1:
            return cls(rate_categories=rate_categories)

        if not rate_model_str.endswith("}"):
            msg = f"Missing end bracket for parameterisation {rate_model_str!r}"
            raise ValueError(msg)

        try:
            alpha = float(rate_model_str[parameter_start + 1 : -1])
        except ValueError:
            msg = f"Parameterisation of Discrete Gamma Model is not a number {rate_model_str!r}"
            raise ValueError(msg) from None

        return cls(rate_categories=rate_categories, alpha=alpha)


class FreeRateModel(RateModel):
    def __init__(self, rate_categories: int | None = None) -> None:
        """FreeRate Model.

        Parameters
        ----------
        rate_categories : int, optional
            The number of rate categories, by default 4.

        References
        ----------
        .. [1] Yang, Ziheng. "A space-time process model for the evolution of
           DNA sequences." Genetics 139.2 (1995): 993-1005.
        .. [2] Soubrier, Julien, et al. "The influence of rate heterogeneity
           among sites on the time dependence of molecular rates." Molecular
           biology and evolution 29.11 (2012): 3345-3358.

        """
        self.rate_categories = rate_categories

    def iqtree_str(self) -> str:
        if self.rate_categories is None:
            return "R"
        return f"R{self.rate_categories}"

    @classmethod
    def from_str(cls, rate_model_str: str) -> "FreeRateModel":
        if len(rate_model_str) == 0:
            msg = "An empty string is not a DiscreteGammaModel."
            raise ValueError(msg)

        if rate_model_str[0] != "R":
            msg = f"A FreeRateModel must start with G but got {rate_model_str!r}."
            raise ValueError(msg)

        if len(rate_model_str) == 1:
            return cls()

        rate_categories = _parse_rate_categories(rate_model_str)
        return cls(rate_categories=rate_categories)


ALL_BASE_RATE_TYPES = [
    RateType(),
    RateType(invariable_sites=True),
    RateType(rate_model=DiscreteGammaModel()),
    RateType(invariable_sites=True, rate_model=DiscreteGammaModel()),
    RateType(rate_model=FreeRateModel()),
    RateType(invariable_sites=True, rate_model=FreeRateModel()),
]

_BASE_RATE_TYPE_DESCRIPTIONS = {
    RateType().iqtree_str(): "no invariable sites, no rate heterogeneity model.",
    RateType(
        invariable_sites=True,
    ).iqtree_str(): "allowing for a proportion of invariable sites.",
    RateType(
        rate_model=DiscreteGammaModel(),
    ).iqtree_str(): "discrete Gamma model (Yang, 1994) with default 4 rate categories. The number of categories can be changed with e.g. +G8.",
    RateType(
        invariable_sites=True,
        rate_model=DiscreteGammaModel(),
    ).iqtree_str(): "invariable site plus discrete Gamma model (Gu et al., 1995).",
    RateType(
        rate_model=FreeRateModel(),
    ).iqtree_str(): "FreeRate model (Yang, 1995; Soubrier et al., 2012) that generalizes the +G model by relaxing the assumption of Gamma-distributed rates. The number of categories can be specified with e.g. +R6 (default 4 categories if not specified). The FreeRate model typically fits data better than the +G model and is recommended for analysis of large data sets.",
    RateType(
        invariable_sites=True,
        rate_model=FreeRateModel(),
    ).iqtree_str(): "invariable site plus FreeRate model.",
}


def get_description(rate_type: RateType) -> str:
    rate_type_str = "".join(c for c in rate_type.iqtree_str() if not c.isdigit())
    return _BASE_RATE_TYPE_DESCRIPTIONS[rate_type_str]


def get_rate_type(
    rate_model: str | RateModel | None = None,
    *,
    invariable_sites: bool | float = False,
) -> RateType:
    """Make a RateType from a chosen rate model and invariable sites.

    Parameters
    ----------
    rate_model : str | RateModel | None, optional
        The chosen rate model, by default None.
    invariable_sites : bool | float, optional
        Invariable sites, by default False.
        If a float in range [0,1) specifies the proportion of invariable sites.

    Returns
    -------
    RateType
        RateType generated from the rate model with invariable sites.

    """
    if isinstance(rate_model, RateModel):
        return RateType(rate_model=rate_model, invariable_sites=invariable_sites)

    if rate_model is None:
        return RateType(invariable_sites=invariable_sites)

    if not isinstance(rate_model, str):
        msg = f"Unexpected type for rate_model: {type(rate_model)}"
        raise TypeError(msg)

    stripped_rate_model = rate_model.lstrip("+")

    rate_model_char = stripped_rate_model[0]
    if rate_model_char == "G":
        rate_model = DiscreteGammaModel.from_str(stripped_rate_model)
    elif rate_model_char == "R":
        rate_model = FreeRateModel.from_str(stripped_rate_model)
    else:
        msg = f"Unexpected value for rate_model {rate_model!r}"
        raise ValueError(msg)

    return RateType(
        rate_model=rate_model,
        invariable_sites=invariable_sites,
    )


def _parse_rate_categories(rate_model_str: str) -> int | None:
    # Assume that the rate model str starts with a G or an R
    parameters_start = rate_model_str.find("{")

    # There is no {} style parameterisation
    if parameters_start == -1:
        categories_str = rate_model_str[1:]
    else:
        categories_str = rate_model_str[1:parameters_start]

    if len(categories_str) == 0:
        return None

    try:
        return int(categories_str)
    except ValueError:
        msg = f"Invalid specification for rate categories {rate_model_str!r}"
        raise ValueError(msg) from None
