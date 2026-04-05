import pandas as pd
import uuid

from immune_nets.models import Repertoire
from immune_nets.creation.algorithms.common_methods import split_tcr_column

class RepertoireFactory:
    @staticmethod
    def fromCSV( name, desc,  path) -> Repertoire:
        df =  pd.read_csv(path)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        new_repertoire = Repertoire(name=name, description=desc)
        new_repertoire.setClones(df)
        new_repertoire.repertoire_id = uuid.uuid4()
        return new_repertoire
