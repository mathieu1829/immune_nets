from src.creation.globalSettings import globalSettings
from src.creation.immuneRepertoire import immuneRepertoire


    #Checks whether creation algorithm is valid
def algorithm(algo):
    def algorithmValidator(strategy=globalSettings().defaultOutputStrategy,**kwargs):
        if type(kwargs['repertoire']) != immuneRepertoire:
            raise TypeError(f"Expected clonotype, got ${type(kwargs['repertoire'])}")
        if not hasattr(kwargs['repertoire'].clones,"name"):
            raise ValueError(f"Clonotype table has no name")
        return strategy(**kwargs, algo = algo)
    return algorithmValidator 

    
