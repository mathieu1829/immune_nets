from abc import ABC, abstractmethod
from immune_nets.entities import GraphStats

class InterGroupStatDistance(ABC):
    @abstractmethod
    def group_dist(self, a:list[GraphStats], b:list[GraphStats]) -> float:
        pass

