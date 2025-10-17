import pandas as pd

from src.models import Repertoire
from src.creation.algorithms.common_methods import split_tcr_column

class RepertoireFactory:
    @staticmethod
    def fromCSV( name, desc,  path) -> Repertoire:
        df =  pd.read_csv(path)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        new_repertoire = Repertoire(name=name, description=desc)
        new_repertoire.setClones(df)
        return new_repertoire
