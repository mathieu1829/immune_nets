import pandas as pd
import uuid

from immune_nets.creation.algorithms.common_methods import split_tcr_column

class ImmuneRepertoire:
    def __init__(self,name,clones,description=" ", repertoire_id=None):
        self.name = name
        self.description = description
        self.clones = clones
        self.repertoire_id = repertoire_id if repertoire_id is not None else uuid.uuid4()


