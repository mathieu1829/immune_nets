import pandas as pd
import uuid

from src.creation.algorithms.common_methods import split_tcr_column

class immuneRepertoire:
    def __init__(self,name,clones,description=" ", repertoire_id=uuid.uuid4):
        self.name = name
        self.description = description
        self.clones = clones
        self.repertoire_id = repertoire_id

    @classmethod
    def fromCSV(cls, name,desc,path):
        df =  pd.read_csv(path)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        new_repertoire = cls(name=name, description=desc, clones=df)
        return new_repertoire
