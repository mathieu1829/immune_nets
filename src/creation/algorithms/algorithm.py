from abc import ABC,abstractmethod
import pandas as pd


    #Checks whether creation algorithm is valid
def algorithm(algo):
    def algorithmValidator(**kwargs):
        strategy = kwargs['strategy']
        assert type(kwargs['clonotypes']) == type(pd.DataFrame()), f"Expected clonotype, got ${type(kwargs['clonotypes'])}"
        assert hasattr(kwargs['clonotypes'],"name") , f"Clonotype table has no name"
        return strategy(**kwargs, algo = algo)
    return algorithmValidator 

    
