from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

def split_tcr_column(x, subunit):
    s1 = x.split(';')
    for sx in s1:
        s = sx.split(':')
        if s[0] == subunit:
            return s[1]

def tcr_alig(x, ref, aligner):
    if x != None:
        alignments = pairwise2.align.localds(x, ref, aligner.substitution_matrix, -10, -1)
        return alignments[0].score
    
