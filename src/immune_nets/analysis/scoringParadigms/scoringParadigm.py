from abc import ABC, abstractmethod

class ScoringParadigm(ABC):
    @abstractmethod
    def compute_score(self, groupedStats) -> float:
        pass
