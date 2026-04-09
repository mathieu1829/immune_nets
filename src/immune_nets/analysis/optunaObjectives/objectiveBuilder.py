from abc import ABC, abstractmethod
from optuna import Trial

class objectiveBuilder(ABC):
    @abstractmethod 
    def __call__(self, trial: Trial) -> float:
        pass
