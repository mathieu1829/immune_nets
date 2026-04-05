from abc import ABC, abstractmethod
from scipy import stats as scipy_stats

class DistributionIdentifier(ABC):
    """
    Interface, base class for distribution identifiers. Enforces implementation of identify_distribution method
    """
    candidates = {
        "powerlaw": scipy_stats.powerlaw,
        "norm": scipy_stats.norm,
        "uniform": scipy_stats.uniform,
        "gamma": scipy_stats.gamma,
    }

    @abstractmethod
    def identify_distribution(self, statDict) -> str:
        """
        Abstract method for identifying provided distribution. 

        :param statDict: dictionairy containing probability mass function, keys are arguments of probability mass function and values of the dictionairy represent values of probability mass function
        :return: name of the distribution
        """

        pass
