from abc import ABC, abstractmethod
from src.entities import GraphStats

class PairwiseStatDistance(ABC):
    @abstractmethod
    def stat_dist(self, a:GraphStats, b:GraphStats) -> float:
        pass

