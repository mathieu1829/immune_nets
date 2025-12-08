from abc import ABC, abstractmethod

class HolisticStatDistance(ABC):
    @abstractmethod
    def compute_holistic_score(self, datasets) -> float:
        pass
