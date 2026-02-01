import pandas as pd

from src.entities import ImmuneRepertoire
from src.creation.algorithms.common_methods import split_tcr_column

class ImmuneRepertoireFactory:
    @staticmethod
    def fromCSV(path, name, desc) -> ImmuneRepertoire:
        df =  pd.read_csv(path)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        new_repertoire = ImmuneRepertoire(name=name, description=desc, clones=df)
        return new_repertoire

    @classmethod
    def fromCSVTest(cls,path) -> ImmuneRepertoire:
        return cls.fromCSV(path=path, name="test repertoire", desc=" ")
