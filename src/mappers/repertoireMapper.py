from src.models import Repertoire
from src.creation.immuneRepertoire import ImmuneRepertoire
import pandas as pd

class RepertoireMapper:
    @staticmethod
    def fromImmuneRepertoire(repertoire: ImmuneRepertoire):
        new_repertoire = Repertoire(repertoire_id=repertoire.repertoire_id,
                             name=repertoire.name,
                             description=repertoire.description
                             )
        new_repertoire.setClones(repertoire.clones)
        return new_repertoire

    @staticmethod
    def toImmuneRepertoire(repertoire: Repertoire):
        clones = pd.DataFrame([{
            "proportion": c.proportion,
            "tcra_aa": c.tcra_aa,
            "tcrb_aa": c.tcrb_aa,
            "cdr3s_nt": c.cdr3s_nt,
            "inkt_evidence": c.inkt_evidence,
            "mait_evidence": c.mait_evidence
        } for c in repertoire.clonotypes])

        return ImmuneRepertoire(repertoire_id=repertoire.repertoire_id,
                                name=repertoire.name,
                                description=repertoire.description,
                                clones=clones
                                )

