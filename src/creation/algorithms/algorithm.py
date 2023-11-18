from src.creation.globalSettings import globalSettings
from src.creation.immuneRepertoire import immuneRepertoire


    #Checks whether creation algorithm is valid
def algorithm(algo):
    def algorithmValidator(strategy=globalSettings().defaultOutputStrategy,**kwargs):
        assert type(kwargs['repertoire']) == immuneRepertoire, f"Expected clonotype, got ${type(kwargs['repertoire'])}"
        assert hasattr(kwargs['repertoire'],"name") , f"Clonotype table has no name"
        return strategy(**kwargs, algo = algo)
    return algorithmValidator 

    
