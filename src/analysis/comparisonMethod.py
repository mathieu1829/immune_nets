from abc import ABC,abstractmethod
import pandas as pd

class comparisonMethod(ABC):

    def __init__(self, strategy, **kwargs):
        self.strategy = strategy
        self.kwargs = kwargs #arguments passed to strategy

    #Checks whether creation algorithm is valid
    def validateComparisonAlgorithm(algo):
        def algorithmValidator(self,**overrideKwargs):
            assert type(overrideKwargs['immuneNet']) == type(pd.DataFrame()), f"Expected DataFrame, got ${type(overrideKwargs['immuneNet'])}"
            assert hasattr(overrideKwargs['immuneNet'],"name") , f"Immune net has no name" 
            mergedKwargs = {**self.kwargs,**overrideKwargs}
            onlyKwargs = {key:self.kwargs[key] for key in self.kwargs if not key in overrideKwargs}
            finalArgs = { key:mergedKwargs[key] for key in mergedKwargs if not key in onlyKwargs }
            return self.strategy(**finalArgs,self = self, algo = algo)
        return algorithmValidator 

    @validateComparisonAlgorithm
    @abstractmethod
    def createGraph(self,**kwargs): #method used to create a graph
        return self.comparisonAlgorithm(**kwargs) #implement this to specify algorithm
    
