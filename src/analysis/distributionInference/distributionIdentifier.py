from abc import ABC, abstractmethod
from scipy import stats as scipy_stats

class DistributionIdentifier(ABC):
    candidates = {
        "powerlaw": scipy_stats.powerlaw,
        "norm": scipy_stats.norm,
        "uniform": scipy_stats.uniform,
        "gamma": scipy_stats.gamma,
    }

    @abstractmethod
    def identify_distribution(self, statDict) -> str:
        pass
