import pandas as pd
import uuid

from src.creation.algorithms.common_methods import split_tcr_column

class ImmuneRepertoire:
    def __init__(self,name,clones,description=" ", repertoire_id=uuid.uuid4):
        self.name = name
        self.description = description
        self.clones = clones
        self.repertoire_id = repertoire_id


