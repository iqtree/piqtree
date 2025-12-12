from abc import ABC, abstractmethod


class IndelDistribution(ABC):
    @abstractmethod
    def iqtree_str(self) -> str:
        pass  # pragma: no cover

    def __str__(self) -> str:
        return self.iqtree_str()


class IndelZipfian(IndelDistribution):
    def __init__(self, exponent: float, max_size: int) -> None:
        self.exponent = exponent
        self.max_size = max_size

    def iqtree_str(self) -> str:
        return f"POW{{{self.exponent}/{self.max_size}}}"


class IndelGeometric(IndelDistribution):
    def __init__(self, mean: float) -> None:
        self.mean = mean

    def iqtree_str(self) -> str:
        return f"GEO{{{self.mean}}}"


class IndelNegativeBinomial(IndelDistribution):
    def __init__(self, mean: float, variance: float) -> None:
        self.mean = mean
        self.variance = variance

    def iqtree_str(self) -> str:
        return f"NB{{{self.mean}/{self.variance}}}"


class IndelLavalette(IndelDistribution):
    def __init__(self, a: float, max_size: int) -> None:
        self.a = a
        self.max_size = max_size

    def iqtree_str(self) -> str:
        return f"LAV{{{self.a}/{self.max_size}}}"
