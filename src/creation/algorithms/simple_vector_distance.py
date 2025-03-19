# return graph (dataframe) based on another df with clonotypes

import numpy as np
from src.creation.algorithms.common_methods import *
from src.creation.algorithms.algorithm import *
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.immuneNetwork import immuneNetwork

@algorithm
def simple_vector_distance( repertoire, distance, threshold = None, **kwargs):
    clonotypes = repertoire.clones
    distanceFun = distance.tcr_dist

    if distance.group == True:
        raise ValueError("Vector distance algorithm can't handle \"group\" distance functions")


    ab_tcr = clonotypes[["tcra_aa", "tcrb_aa"]].dropna();
    tcr_npa = clonotypes[["tcra_aa", "tcrb_aa"]].dropna().to_numpy()
    dist_al_trcb = np.zeros(np.shape(tcr_npa)[0] * np.shape(tcr_npa)[0]).reshape(np.shape(tcr_npa)[0],
                                                                                 np.shape(tcr_npa)[0])
    dist_al_trcb += -1

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
    
    ab_tcr['alpha_closest_to_consensus'] = ab_tcr['tcra_aa'].apply(lambda x: distanceFun(x,consensus_alpha_seq))
    ab_tcr['beta_closest_to_consensus'] = ab_tcr['tcrb_aa'].apply(lambda x: distanceFun(x,consensus_beta_seq))


    closest_alpha = ab_tcr[['tcra_aa','alpha_closest_to_consensus']].sort_values('alpha_closest_to_consensus').dropna().to_numpy()[:6,0]
    closest_beta = ab_tcr[['tcrb_aa','beta_closest_to_consensus']].sort_values('beta_closest_to_consensus').dropna().to_numpy()[:6,0]

    for idx,amino in enumerate(closest_alpha):
        ab_tcr[f'a{idx}'] = ab_tcr['tcra_aa'].apply(lambda x : distanceFun(x,amino))

    for idx,amino in enumerate(closest_beta):
        ab_tcr[f'b{idx}'] = ab_tcr['tcrb_aa'].apply(lambda x : distanceFun(x,amino))

    dist_mat = pd.DataFrame(
    squareform(pdist(ab_tcr[["a0","b0", "a1", "b1", "a2", "b2", "a3", "b3", "a4", "b4", "a5", "b5"]])),
    columns = ab_tcr.index,
    index = ab_tcr.index
    ) 
    
    threshold = np.nanmean(dist_mat.to_numpy()) / 4 if threshold is None else threshold
    # threshold = np.mean(dist_mat.fillna(0).to_numpy()) / 4 if threshold is None else threshold #produces smaller threshold
    dist_mat = np.tril(dist_mat, k=-1)
    # threshold = np.nanmean(dist_mat) / 4 if threshold is None else threshold #produces even smaller threshold
    for i in range(len(dist_mat)):
        for j in range(i,len(dist_mat)):
            dist_mat[i,j] = float('INF')
    # print(dist_mat)
    # print(f"threshold: {threshold}")
    matrix_cutoff = np.where(dist_mat < threshold)

    d = {'r1': matrix_cutoff[0], 'r2': matrix_cutoff[1]}
    df_net = pd.DataFrame(data=d)
    df_net.name = clonotypes.name
    immuneNet = immuneNetwork(df_net, "simple_vector_distance", np.unique(repertoire.clones["sampleID"].to_numpy()),str(distance) , threshold , len(ab_tcr) ) 

    return immuneNet


