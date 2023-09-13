from abc import ABC,abstractmethod
import pandas as pd

class algorithm(ABC):

    def __init__(self, strategy, **kwargs):
        self.strategy = strategy
        self.kwargs = kwargs #arguments passed to strategy

    #Checks whether creation algorithm is valid
    def validateCreationAlgorithm(algo):
        def algorithmValidator(self,**overrideKwargs):
            assert type(overrideKwargs['clonotypes']) == type(pd.DataFrame()), f"Expected clonotype, got ${type(overrideKwargs['clonotypes'])}"
            assert hasattr(overrideKwargs['clonotypes'],"name") , f"Clonotype table has no name"
            mergedKwargs = {**self.kwargs,**overrideKwargs}
            onlyKwargs = {key:self.kwargs[key] for key in self.kwargs if not key in overrideKwargs}
            finalArgs = { key:mergedKwargs[key] for key in mergedKwargs if not key in onlyKwargs }
            return self.strategy(**finalArgs,self = self, algo = algo)
        return algorithmValidator 

    @validateCreationAlgorithm
    @abstractmethod
    def createGraph(self,**kwargs): #method used to create a graph
        return self.creationAlgorithm(**kwargs) #implement this to specify algorithm
    
