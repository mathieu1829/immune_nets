from abc import ABC, abstractmethod
from optuna import Trial

class objectiveBuilder(ABC):
    '''
    Interface for building Optuna objective functions.

    Enables defining objective functions as callable objects, allowing
    additional parameters to be configured via the constructor while still
    exposing a ``__call__`` method compatible with Optuna.
    '''
    @abstractmethod 
    def __call__(self, trial: Trial) -> float:
        """
        Evaluate the objective function for a given Optuna trial.

        :param trial: Optuna trial object used to suggest parameters.
        :return: Objective value to be optimized.
        """
        pass
