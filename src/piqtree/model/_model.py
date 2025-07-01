from piqtree.model._freq_type import FreqType, get_freq_type
from piqtree.model._rate_type import RateModel, get_rate_type
from piqtree.model._substitution_model import SubstitutionModel, get_substitution_model


class Model:
    """Specification for substitution models.

    Stores the substitution model with base frequency settings.
    """

    def __init__(
        self,
        submod_type: str | SubstitutionModel,
        freq_type: str | FreqType | None = None,
        rate_model: str | RateModel | None = None,
        *,
        invariable_sites: bool | float = False,
    ) -> None:
        """Construct Model class.

        Parameters
        ----------
        submod_type : str | SubstitutionModel
            The substitution model to use
        freq_type : str | FreqType | None, optional
            State frequency specification, by default None. (defaults
            to empirical base frequencies if not specified by model).
        rate_model : str | RateModel | None, optional
            Rate heterogeneity across sites model, by default
            no Gamma, and no FreeRate.
        invariable_sites : bool | float, optional
            Invariable sites, by default False.
            If a float in range [0,1) specifies the proportion of invariable sites.

        """
        self.submod_type = get_substitution_model(submod_type)
        self.freq_type = get_freq_type(freq_type) if freq_type else None
        self.rate_type = (
            get_rate_type(rate_model, invariable_sites=invariable_sites)
            if rate_model is not None or invariable_sites
            else None
        )

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        attrs = [
            f"submod_type={self.submod_type.iqtree_str()}",
            f"freq_type={getattr(self.freq_type, 'name', None)}",
            f"rate_type={getattr(self.rate_type, 'name', None)}",
        ]
        return f"Model({', '.join(attrs)})"

    def __str__(self) -> str:
        """Convert the model into the IQ-TREE representation.

        Returns
        -------
        str
            The IQ-TREE representation of the mode.

        """
        iqtree_extra_args = (
            x for x in (self.freq_type, self.rate_type) if x is not None
        )
        return "+".join(x.iqtree_str() for x in [self.submod_type, *iqtree_extra_args])

    @property
    def rate_model(self) -> RateModel | None:
        """The RateModel used, if one is chosen.

        Returns
        -------
        RateModel | None
            The RateModel used by the Model.

        """
        return self.rate_type.rate_model if self.rate_type else None

    @property
    def invariable_sites(self) -> bool:
        """Whether invariable sites are used.

        Returns
        -------
        bool
            True if invariable sites are used by the model, False otherwise.

        """
        return self.rate_type.invariable_sites if self.rate_type else False

    @property
    def proportion_invariable_sites(self) -> float | None:
        """The proportion of invariable sites if specified.

        Returns
        -------
        float | None
            The proportion of invariable sites if specified, None otherwise.

        """
        return self.rate_type.proportion_invariable if self.rate_type else None


def make_model(iqtree_str: str) -> Model:
    """Convert an IQ-TREE model specification into a Model class.

    Parameters
    ----------
    iqtree_str : str
        The IQ-TREE model string.

    Returns
    -------
    Model
        The equivalent Model class.
    """
    if "+" not in iqtree_str:
        return Model(iqtree_str)

    sub_mod_str, components = iqtree_str.split("+", maxsplit=1)

    freq_type = None
    invariable_sites: float | bool | None = None
    rate_model = None

    for component in components.split("+"):
        if component.startswith("F"):
            if freq_type is not None:
                msg = f"Model {iqtree_str!r} contains multiple base frequency specifications."
                raise ValueError(msg)
            freq_type = component
        elif component.startswith("I"):
            if invariable_sites is not None:
                msg = f"Model {iqtree_str!r} contains multiple specifications for invariable sites."
                raise ValueError(msg)
            invariable_sites = _parse_invariable_sites(component)

        elif component.startswith(("G", "R")):
            if rate_model is not None:
                msg = f"Model {iqtree_str!r} contains multiple rate heterogeneity specifications."
                raise ValueError(msg)
            rate_model = component
        else:
            msg = f"Model {iqtree_str!r} contains unexpected component."
            raise ValueError(msg)

    if invariable_sites is None:
        invariable_sites = False

    return Model(sub_mod_str, freq_type, rate_model, invariable_sites=invariable_sites)


def _parse_invariable_sites(component: str) -> bool | float:
    # Assumes that component starts with "I"
    remainder = component[1:]
    if len(remainder) == 0:
        return True

    if remainder[0] != "{" or remainder[-1] != "}":
        msg = f"Invalid specification for proportion of invariable sites, got '{component}'."
        raise ValueError(msg)

    number_part = remainder[1:-1]
    try:
        return float(number_part)
    except ValueError:
        msg = f"Failed to read proportion of invariable sites, got '{component}'"
        raise ValueError(msg) from None
