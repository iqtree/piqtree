# Use different kinds of substitution models

piqtree currently supports all named IQ-TREE DNA models including Lie Markov models, empirical amino-acid exchange rate matrices,
as well as specification for base frequencies and rate heterogeneity across sites.

We use the [`Model`](../api/model/Model.md) class to represent the substitution model which can be constructed from strings, or
using enums. Most functions can already compute this directly from the passed model string, so usage of [`make_model`](../api/model/make_model.md) may be preferable. Substitution models can be combined with specification for base frequencies, and rate heterogeneity across sites models.

## Usage

### Standard DNA Models

DNA models may be specified using the [`StandardDnaModel`](../api/model/SubstitutionModel.md#piqtree.model.StandardDnaModel) enum, or by using
the IQ-TREE string representation. A full list of supported DNA models is available [here](../api/model/SubstitutionModel.md#piqtree.model.StandardDnaModel).

```python
from piqtree import Model
from piqtree.model import StandardDnaModel

hky_model_1 = Model("HKY")
hky_model_2 = Model(StandardDnaModel.HKY)
```

#### Parameterisation of Standard DNA Models

IQ-TREE allows DNA models to be parameterised by specifying relative substitution rates. See the [IQ-TREE documentation](https://iqtree.github.io/doc/Substitution-Models) for the precise definition of the parameters for the chosen model.

```python
from piqtree import Model
from piqtree.model import StandardDnaModel

model = Model("GTR{1.0,2.0,1.5,3.7,2.8}")
model = Model(StandardDnaModel.GTR([1.0, 2.0, 1.5, 3.7, 2.8]))
```

### Lie Markov Models

Lie Markov models may be specified using the [`LieModel`](../api/model/SubstitutionModel.md#piqtree.model.LieModel) enum, or by using
the IQ-TREE string representation. A full list of supported DNA models is available [here](../api/model/SubstitutionModel.md#piqtree.model.LieModel). The pairing prefix may also be specified - RY for purine-pyrimidine pairing (default); WS for weak-strong pairing; and MK aMino-Keto pairing.

```python
from piqtree import Model
from piqtree.model import LieModel

lie_ws_6_6_model_1 = Model("WS6.6")
lie_ws_6_6_model_2 = Model(LieModel.LIE_6_6("WS"))

lie_12_12_model_1 = Model("12.12")
lie_12_12_model_2 = Model(LieModel.LIE_12_12)
```

#### Parameterisation of Lie Markov Models

IQ-TREE supports the parameterisation of Lie Markov models. The first number (before the period) corresponds to the basis matrices used by the model. The first basis matrix is unparameterised, and the parameters refer to the weights (bounded between -0.98 and 0.98 exclusive) of the remaining basis matrices of the model.

```python
from piqtree import Model
from piqtree.model import LieModel

model = Model("MK3.3b{0.3,-0.5}")
model = Model(LieModel.LIE_3_3b("MK", [0.3, -0.5]))
```

### Amino-acid Models

Amino-acid models may be specified using the [`AaModel`](../api/model/SubstitutionModel.md#piqtree.model.AaModel) enum, or by using
the IQ-TREE string representation. A full list of supported amino-acid models is available [here](../api/model/SubstitutionModel.md#piqtree.model.AaModel).

```python
from piqtree import Model
from piqtree.model import AaModel

dayhoff_model_1 = Model("Dayhoff")
dayhoff_model_2 = Model(AaModel.Dayhoff)

nq_yeast_model_1 = Model("NQ.yeast")
nq_yeast_model_2 = Model(AaModel.NQ_yeast)
```

### Base Frequencies

Three types of base frequencies can be specified using the [`FreqType`](../api/model/FreqType.md) enum, or additionally the [`CustomBaseFreq`](../api/model/FreqType.md) for fixing base frequencies. Alternatively the IQ-TREE string representation can be used.

If not specified, the chosen model's default settings are used.

- [`F`](../api/model/FreqType.md#piqtree.model.FreqType.F): Empirical base frequencies. String representation also used to fix base frequencies.
- [`FQ`](../api/model/FreqType.md#piqtree.model.FreqType.FQ): Equal base frequencies.
- [`FO`](../api/model/FreqType.md#piqtree.model.FreqType.FO): Optimised base frequencies by maximum-likelihood.

```python
from piqtree import Model
from piqtree.model import CustomBaseFreq, FreqType

# Default for the GTR model
empirical_freqs_1 = Model("GTR", freq_type="F")
empirical_freqs_2 = Model("GTR", freq_type=FreqType.F)

equal_freqs_1 = Model("GTR", freq_type="FQ")
equal_freqs_2 = Model("GTR", freq_type=FreqType.FQ)

opt_freqs_1 = Model("GTR", freq_type="FO")
opt_freqs_2 = Model("GTR", freq_type=FreqType.FO)

custom_freqs_1 = Model("GTR", freq_type="F{0.1,0.2,0.3,0.4}")
custom_freqs_2 = Model("GTR", freq_type=CustomBaseFreq([0.1, 0.2, 0.3, 0.4]))
```

### Rate Heterogeneity

#### Invariable Sites

A boolean flag can be specified when constructing the [`Model`](../api/model/Model.md) class to allow for a proportion of invariable sites. To fix the proportion of invariable sites, a float representing the proportion can be used instead.

```python
from piqtree import Model

without_invar_sites = Model("TIM", invariable_sites=False) # Default
with_invar_sites = Model("TIM", invariable_sites=True)

prop_invar = Model("TIM", invariable_sites=0.1)
```

#### Discrete Gamma Model

We support the [`DiscreteGammaModel`](../api/model/RateModel.md#piqtree.model.DiscreteGammaModel) allowing for a variable number of rate categories (by default 4). The Gamma shape parameter (alpha) can also be set.

```python
from piqtree import Model
from piqtree.model import DiscreteGammaModel

# 4 rate categories, no invariable sites
k81_discrete_gamma_4 = Model("K81", rate_model=DiscreteGammaModel())

# 8 rate categories, with invariable sites
k81_invar_discrete_gamma_8 = Model("K81", rate_model=DiscreteGammaModel(8), invariable_sites=True)

# 4 rate categories, alpha=0.2
k81_alpha = Model("K81", rate_model=DiscreteGammaModel(alpha=0.2))
```

#### FreeRate Model

We support the [`FreeRateModel`](../api/model/RateModel.md#piqtree.model.FreeRateModel) allowing for a variable number of rate categories (by default 4). The FreeRate weights for each of the rate categories and the corresponding rates may also be fixed.

```python
from piqtree import Model
from piqtree.model import FreeRateModel

# 4 rate categories, no invariable sites
sym_free_rate_4 = Model("SYM", rate_model=FreeRateModel())

# 8 rate categories, with invariable sites
sym_invar_free_rate_8 = Model("SYM", rate_model=FreeRateModel(8), invariable_sites=True)

# 2 rate categories with specified weights and rates
sym_free_rate_parameterised = Model("SYM", rate_model=FreeRateModel(2, weights=[0.2, 0.8], rates=[2.5, 0.625]))
```

### Making Model Classes from IQ-TREE Strings

For the supported model types, the Model class can be created by using [`make_model`](../api/model/make_model.md) on the IQ-TREE string representation of the model.

```python
from piqtree import make_model

model = make_model("GTR+FQ+I+R3")
```

## See also

- Use a [`Model`](../api/model/Model.md) to construct a maximum likelihood tree: ["Construct a maximum likelihood phylogenetic tree"](construct_ml_tree.md).
- Use a [`Model`](../api/model/Model.md) to fit branch lengths to a tree topology: ["Fit branch lengths to a tree topology from an alignment"](fit_tree_topology.md).
