from abc import ABC,abstractmethod
import pandas as pd


 
 
class algorithm(ABC):

    def __init__(self, strategy, **kwargs):
        self.strategy = strategy
        self.kwargs = kwargs #arguments passed to strategy

    #Checks whether creation algorithm is valid
    def validateCreationAlgorithm(algo):
        def algorithmValidator(self,clonotypes,**overrideKwargs):
            assert type(clonotypes) == type(pd.DataFrame()), f"Expected clonotype, got ${type(clonotypes)}"
            assert hasattr(clonotypes,"name") , f"Clonotype table has no name"
            return self.strategy(**(overrideKwargs if bool(overrideKwargs) else self.kwargs),self = self, algo = algo,clonotypes=clonotypes)
        return algorithmValidator 

    @validateCreationAlgorithm
    @abstractmethod
    def createGraph(self,clonotypes,**kwargs): #method used to create a graph
        return self.creationAlgorithm(clonotypes, **kwargs) #implement this to specify algorithm
    
