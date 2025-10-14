from src.creation.globalSettings import globalSettings

from src.creation.immuneRepertoire import ImmuneRepertoire

    #Checks whether creation algorithm is valid
def algorithm(algo):
    def algorithmValidator(**kwargs):
        if type(kwargs['repertoire']) != ImmuneRepertoire:
            raise TypeError(f"Expected clonotype, got ${type(kwargs['repertoire'])}")
        return algo(**kwargs)
    return algorithmValidator 

    
