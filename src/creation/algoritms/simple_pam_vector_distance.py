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

    unique_amino_acids = np.array(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V',''])
    
    

    alpha_profile = []
    beta_profile = []

    for seq in range(tcr_npa.shape[0]):
        alpha_seq = list(tcr_npa[seq,0])
        beta_seq = list(tcr_npa[seq,1])
        for idx,amino in enumerate(alpha_seq):
            if idx == len(alpha_profile):
                alpha_profile.append(np.zeros(unique_amino_acids.shape[0]))
                alpha_profile[-1][-1] += seq
            alpha_profile[idx][np.where(unique_amino_acids == amino)]+=1
        for i in range(len(alpha_profile) - (len(alpha_seq))):
            alpha_profile[len(alpha_seq) - 1 + i][-1]+=1
        for idx,amino in enumerate(beta_seq):
            if idx == len(beta_profile):
                beta_profile.append(np.zeros(unique_amino_acids.shape[0]))
                beta_profile[-1][-1] += seq
            beta_profile[idx][np.where(unique_amino_acids == amino)]+=1
        for i in range(len(beta_profile) - (len(beta_seq))):
            beta_profile[len(beta_seq) - 1 + i][-1]+=1



    consensus_alpha_arr = np.array([unique_amino_acids[alpha_profile[peptide].argmax()] for peptide in range(len(alpha_profile))])
    consensus_beta_arr = np.array([unique_amino_acids[beta_profile[peptide].argmax()] for peptide in range(len(beta_profile))])

    consensus_alpha_seq = ""
    for i in range(len(alpha_profile)):
        consensus_alpha_seq+=consensus_alpha_arr[i]

    consensus_beta_seq = ""
    for i in range(len(beta_profile)):
        consensus_beta_seq+=consensus_beta_arr[i]
    
    df['alpha_closest_to_consensus'] = df['tcra_aa'].apply(lambda x: tcr_alig(x,consensus_alpha_seq))
    df['beta_closest_to_consensus'] = df['tcrb_aa'].apply(lambda x: tcr_alig(x,consensus_beta_seq))


    closest_alpha = df[['tcra_aa','alpha_closest_to_consensus']].sort_values('alpha_closest_to_consensus').dropna().to_numpy()[-6:,0]
    closest_beta = df[['tcrb_aa','beta_closest_to_consensus']].sort_values('beta_closest_to_consensus').dropna().to_numpy()[-6:,0]

    for idx,amino in enumerate(closest_alpha):
        df[f'a{idx}'] = df['tcra_aa'].apply(lambda x : tcr_alig(x,amino))

    for idx,amino in enumerate(closest_beta):
        df[f'b{idx}'] = df['tcrb_aa'].apply(lambda x : tcr_alig(x,amino))

    dist_mat = pd.DataFrame(
    squareform(pdist(df[["a0","b0", "a1", "b1", "a2", "b2", "a3", "b3", "a4", "b4", "a5", "b5"]])),
    columns = df.index,
    index = df.index
    )    


    treshold = np.mean(dist_mat) / 4# todo -> treat as parameter
    dist_mat = np.tril(dist_mat, k=-1)
    matrix_cutoff = np.where(dist_mat > treshold)

    d = {'r1': matrix_cutoff[0], 'r2': matrix_cutoff[1]}
    df_net = pd.DataFrame(data=d)
    return df_net


