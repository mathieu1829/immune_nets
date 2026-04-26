from abc import ABC, abstractmethod
from immune_nets.entities.graphStats import GraphStats

class ScoringParadigm(ABC):
    """
    Interface for scoring paradigms. Enforces implementation of the compute_score method.
    """
    @abstractmethod
    def compute_score(self, groupedStats: dict[str, list[GraphStats]]) -> float:
        """
        Compute a distance (or score) between groups of GraphStats.

        :param grouped_stats: A dictionary mapping disease state names to lists of
                              GraphStats objects representing networks build on samples
                              derived from patients in that disease state.
        :return: A float representing the computed score
        """
        pass
