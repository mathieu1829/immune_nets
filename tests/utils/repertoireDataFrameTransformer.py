import pandas as pd
from tests.utils.repertoireData import repertoireData
import numpy as np
from Bio.Seq import Seq
from Bio.Data.CodonTable import standard_dna_table

class repertoireDataFrameTransformer:
    def reverse_translate_simple(self,protein_seq):
        aa_to_codon = {}
        for codon, aa in standard_dna_table.forward_table.items():
            if aa not in aa_to_codon:
                aa_to_codon[aa] = codon  # Pick the first codon for each amino acid
        aa_to_codon['*'] = "TAA"  # Stop codon

        # Translate the protein sequence
        try:
            return "".join(aa_to_codon[aa] for aa in protein_seq)
        except KeyError as e:
            raise ValueError(f"Invalid amino acid in protein sequence: {e.args[0]}") from e

    def transform(self, repertoire: repertoireData) -> pd.DataFrame:
        data = {
                "clonotype_id": [ f"clonotype{i}" for i in range(len(repertoire.tcrA))],
                "frequency": repertoire.frequency,
                "proportion": np.array([ i/repertoire.frequency.sum() for i in repertoire.frequency]),
                "cdr3s_aa": np.array([ f"TRB:{tcrB};TRA:{tcrA}" for tcrA,tcrB in zip(repertoire.tcrA,repertoire.tcrB)]),
                "cdr3s_nt": np.array([ f"TRB:{self.reverse_translate_simple(tcrB)};TRA:{self.reverse_translate_simple(tcrA)}" for tcrA,tcrB in zip(repertoire.tcrA,repertoire.tcrB)]),
                "inkt_evidence": [ 1 for i in range(len(repertoire.tcrA))], 
                "mait_evidence": [ 1 for i in range(len(repertoire.tcrA))] 
                } 
        df = pd.DataFrame(data=data)
        return df
