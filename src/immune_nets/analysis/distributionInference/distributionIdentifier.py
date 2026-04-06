from abc import ABC, abstractmethod
from scipy import stats as scipy_stats

class DistributionIdentifier(ABC):
    """
    Base class for distribution identifiers.

    Provides an interface for identifying a probability distribution based on
    input data, and defines a set of candidate distributions.
    """
    candidates = {
        "powerlaw": scipy_stats.powerlaw,
        "normal": scipy_stats.norm,
        "uniform": scipy_stats.uniform,
        "gamma": scipy_stats.gamma,
    }
    """Mapping of distribution names to corresponding ``scipy.stats`` objects.

    Defines the set of candidate distributions used during identification.
    """

    @abstractmethod
    def identify_distribution(self, statDict) -> str:
        """
        Identify the distribution represented by the given data.

        :param statDict: Dictionary representing a probability mass function,
            where keys are input values and values are the corresponding
            probabilities.
        :return: Name of the identified distribution.
        """

        pass
