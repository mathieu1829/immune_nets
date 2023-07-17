# return graph (dataframe) based on another df with clonotypes

import numpy as np
from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from common_methods import *
import logging
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("PAM250")


def tcr_alig(x, ref):
    if x != None:
        alignments = pairwise2.align.localds(x, ref, aligner.substitution_matrix, -10, -1)
        return alignments[0].score


def create_PAM250_vector_network(df):
    tcr_npa = df[["tcra_aa", "tcrb_aa"]].dropna().to_numpy()
    dist_al_trcb = np.zeros(np.shape(tcr_npa)[0] * np.shape(tcr_npa)[0]).reshape(np.shape(tcr_npa)[0],
                                                                                 np.shape(tcr_npa)[0])
    dist_al_trcb += -1
    matrix_len = len(dist_al_trcb)
    longest_alpha_seq = np.array([len(tcr_npa[i,0]) for i in range(tcr_npa.shape[0])]).max()
    longest_beta_seq = np.array([len(tcr_npa[i,1]) for i in range(tcr_npa.shape[0])]).max()

    alpha_mat = np.array([ [ tcr_npa[seq,0][peptide] if peptide < len(tcr_npa[seq,0]) else '' for peptide in range(longest_alpha_seq)] for seq in range(tcr_npa.shape[0])])
    beta_mat = np.array([ [ tcr_npa[seq,1][peptide] if peptide < len(tcr_npa[seq,1]) else '' for peptide in range(longest_beta_seq)] for seq in range(tcr_npa.shape[1])])

    alpha_unique = np.unique((alpha_mat).tolist())
    beta_unique = np.unique((beta_mat).tolist())

    alpha_profile = np.array([[(alpha_mat[:,peptide_num] == peptide_symbol).sum()  for peptide_num in range(longest_alpha_seq)] for peptide_symbol in alpha_unique])
    beta_profile = np.array([[(beta_mat[:,peptide_num] == peptide_symbol).sum()  for peptide_num in range(longest_beta_seq)] for peptide_symbol in beta_unique])

    consensus_alpha_arr = np.array([alpha_unique[(alpha_profile[:,peptide]).argmax()] for peptide in range(longest_alpha_seq)])
    consensus_beta_arr = np.array([beta_unique[(beta_profile[:,peptide]).argmax()] for peptide in range(longest_beta_seq)])

    consensus_alpha_seq = ""
    for i in range(longest_alpha_seq):
        consensus_alpha_seq+=consensus_alpha_arr[i]
    
    consensus_beta_seq = ""
    for i in range(longest_beta_seq):
        consensus_beta_seq+=consensus_beta_arr[i]

    df['alpha_closest_to_consensus'] = df['tcra_aa'].apply(lambda x: tcr_alig(x,consensus_alpha_seq))
    df['beta_closest_to_consensus'] = df['tcrb_aa'].apply(lambda x: tcr_alig(x,consensus_beta_seq))

    closest_alpha_idx = np.array([i for i in np.sort(df['alpha_closest_to_consensus'].dropna().to_numpy())[-6:]])
    # closest_beta_idx = np.array([np.where(df['beta_closest_to_consensus'].to_numpy() == i)[0] for i in np.sort(df['beta_closest_to_consensus'].dropna().to_numpy())[-6:]])
    print(closest_alpha_idx)

    for idx,val in enumerate(closest_alpha_idx):
        # df[f'a{idx}'] = df['tcra_aa'].apply(lambda x : tcr_alig(x,df['tcra_aa'][val]))
        print(f"{idx} {val}")

    # for idx,val in enumerate(closest_beta_idx):
    #     df[f'b{idx}'] = df['tcra_aa'].apply(lambda x : tcr_alig(x,df['tcrb_aa'][val]))
    #
    # dist_mat = pd.DataFrame(
    # squareform(pdist(df[["a0","b0", "a1", "b1", "a2", "b2", "a3", "b3", "a4", "b4", "a5", "b5"]])),
    # columns = df.index,
    # index = df.index
    # )    


df = pd.read_csv("../../../tests/test_data/test_clonotypes.csv")
df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))




create_PAM250_vector_network(df)

