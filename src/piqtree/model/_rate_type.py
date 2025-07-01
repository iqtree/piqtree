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
        invariant_sites: bool | float = False,
        rate_model: RateModel | None = None,
    ) -> None:
        """Rate heterogeneity across sites model.

        Parameters
        ----------
        invariant_sites : bool | float, optional
            Invariable sites, by default False.
            If a float in range [0,1) specifies the proportion of invariable sites.
        rate_model : RateModel | None, optional
            Discrete Gamma Model or FreeRate Model, by default None.

        """

        if not isinstance(invariant_sites, bool):
            if not (0 <= invariant_sites < 1):
                msg = "The proportion of invaraint sites must be in the range [0,1)"
                raise ValueError(msg)
            self.invariant_sites = True
            self.proportion_invariant: float | None = invariant_sites
        else:
            self.invariant_sites = invariant_sites
            self.proportion_invariant = None

        self.rate_model = rate_model

    def iqtree_str(self) -> str:
        """Convert to an iqtree settings string.

        Returns
        -------
        str
            String parsable by IQ-TREE for the rate heterogeneity model.

        """
        rate_type_str = "I" if self.invariant_sites else ""
        if self.proportion_invariant is not None:
            rate_type_str += f"{{{self.proportion_invariant}}}"

        if self.rate_model is None:
            return rate_type_str
        # Invariant sites and model need to be joined by a '+'
        if self.invariant_sites:
            rate_type_str += "+"
        return rate_type_str + self.rate_model.iqtree_str()

    @property
    def name(self) -> str:
        return self.iqtree_str()


class DiscreteGammaModel(RateModel):
    def __init__(self, rate_categories: int | None = None) -> None:
        """Discrete Gamma Model.

        Parameters
        ----------
        rate_categories : int, optional
            The number of rate categories, by default 4.

        References
        ----------
        .. [1] Yang, Ziheng. "Maximum likelihood phylogenetic estimation from
           DNA sequences with variable rates over sites: approximate methods."
           Journal of Molecular evolution 39 (1994): 306-314.

        """
        self.rate_categories = rate_categories

    def iqtree_str(self) -> str:
        if self.rate_categories is None:
            return "G"
        return f"G{self.rate_categories}"


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


ALL_BASE_RATE_TYPES = [
    RateType(),
    RateType(invariant_sites=True),
    RateType(rate_model=DiscreteGammaModel()),
    RateType(invariant_sites=True, rate_model=DiscreteGammaModel()),
    RateType(rate_model=FreeRateModel()),
    RateType(invariant_sites=True, rate_model=FreeRateModel()),
]

_BASE_RATE_TYPE_DESCRIPTIONS = {
    RateType().iqtree_str(): "no invariable sites, no rate heterogeneity model.",
    RateType(
        invariant_sites=True,
    ).iqtree_str(): "allowing for a proportion of invariable sites.",
    RateType(
        rate_model=DiscreteGammaModel(),
    ).iqtree_str(): "discrete Gamma model (Yang, 1994) with default 4 rate categories. The number of categories can be changed with e.g. +G8.",
    RateType(
        invariant_sites=True,
        rate_model=DiscreteGammaModel(),
    ).iqtree_str(): "invariable site plus discrete Gamma model (Gu et al., 1995).",
    RateType(
        rate_model=FreeRateModel(),
    ).iqtree_str(): "FreeRate model (Yang, 1995; Soubrier et al., 2012) that generalizes the +G model by relaxing the assumption of Gamma-distributed rates. The number of categories can be specified with e.g. +R6 (default 4 categories if not specified). The FreeRate model typically fits data better than the +G model and is recommended for analysis of large data sets.",
    RateType(
        invariant_sites=True,
        rate_model=FreeRateModel(),
    ).iqtree_str(): "invariable site plus FreeRate model.",
}


def get_description(rate_type: RateType) -> str:
    rate_type_str = "".join(c for c in rate_type.iqtree_str() if not c.isdigit())
    return _BASE_RATE_TYPE_DESCRIPTIONS[rate_type_str]


def get_rate_type(
    rate_model: str | RateModel | None = None,
    *,
    invariant_sites: bool | float = False,
) -> RateType:
    """Make a RateType from a chosen rate model and invariant sites.

    Parameters
    ----------
    rate_model : str | RateModel | None, optional
        The chosen rate model, by default None.
    invariant_sites : bool | float, optional
        Invariable sites, by default False.
        If a float in range [0,1) specifies the proportion of invariable sites.

    Returns
    -------
    RateType
        RateType generated from the rate model with invariant sites.

    """
    if isinstance(rate_model, RateModel):
        return RateType(rate_model=rate_model, invariant_sites=invariant_sites)

    if rate_model is None:
        return RateType(invariant_sites=invariant_sites)

    if not isinstance(rate_model, str):
        msg = f"Unexpected type for rate_model: {type(rate_model)}"
        raise TypeError(msg)

    stripped_rate_model = rate_model.lstrip("+")
    if len(stripped_rate_model) == 1:
        rate_categories = None
    else:
        integer_part = stripped_rate_model[1:]
        if not integer_part.isdigit():
            msg = f"Unexpected value for rate_model {rate_model!r}"
            raise ValueError(msg)

        rate_categories = int(integer_part)

    if stripped_rate_model[0] == "G":
        return RateType(
            rate_model=DiscreteGammaModel(rate_categories=rate_categories),
            invariant_sites=invariant_sites,
        )

    if stripped_rate_model[0] == "R":
        return RateType(
            rate_model=FreeRateModel(rate_categories=rate_categories),
            invariant_sites=invariant_sites,
        )

    msg = f"Unexpected value for rate_model {rate_model!r}"
    raise ValueError(msg)
