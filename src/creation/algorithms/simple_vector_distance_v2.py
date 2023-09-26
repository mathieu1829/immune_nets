# return graph (dataframe) based on another df with clonotypes
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.algorithms.common_methods import *
from src.creation.algorithms.algorithm import *

class simple_vector_distance_v2(algorithm):
    aligner = PairwiseAligner()
    def creationAlgorithm(self,clonotypes, matrix, threshold_alfa = None, threshold_beta = None, **kwargs):
        self.aligner.substitution_matrix = substitution_matrices.load(matrix)
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
        
        clonotypes['alpha_closest_to_consensus'] = clonotypes['tcra_aa'].apply(lambda x: tcr_alig(x,consensus_alpha_seq,self.aligner))
        clonotypes['beta_closest_to_consensus'] = clonotypes['tcrb_aa'].apply(lambda x: tcr_alig(x,consensus_beta_seq,self.aligner))


        closest_alpha = clonotypes[['tcra_aa','alpha_closest_to_consensus']].sort_values('alpha_closest_to_consensus').dropna().to_numpy()[-6:,0]
        closest_beta = clonotypes[['tcrb_aa','beta_closest_to_consensus']].sort_values('beta_closest_to_consensus').dropna().to_numpy()[-6:,0]

        for idx,amino in enumerate(closest_alpha):
            clonotypes[f'a{idx}'] = clonotypes['tcra_aa'].apply(lambda x : tcr_alig(x,amino,self.aligner))

        for idx,amino in enumerate(closest_beta):
            clonotypes[f'b{idx}'] = clonotypes['tcrb_aa'].apply(lambda x : tcr_alig(x,amino,self.aligner))

        dist_mat_alpha = pd.DataFrame(
        squareform(pdist(clonotypes[["a0", "a1", "a2", "a3", "a4", "a5"]])),
        columns = clonotypes.index,
        index = clonotypes.index
        ) 
        
        dist_mat_beta = pd.DataFrame(
        squareform(pdist(clonotypes[["b0", "b1", "b2", "b3", "b4", "b5"]])),
        columns = clonotypes.index,
        index = clonotypes.index
        ) 


        treshold_alpha = np.nanmean(dist_mat_alpha.to_numpy()) / 4 if threshold_alfa is None else threshold_alfa# todo -> treat as parameter, and in this case we could make different tresholds for alpha and beta
        treshold_beta = np.nanmean(dist_mat_beta.to_numpy())  / 4 if threshold_beta is None else threshold_beta# todo -> treat as parameter, and in this case we could make different tresholds for alpha and beta
        dist_mat_alpha = np.tril(dist_mat_alpha, k=-1)
        dist_mat_beta = np.tril(dist_mat_beta, k=-1)
        matrix_cutoff = np.where(((dist_mat_alpha < treshold_alpha) & (dist_mat_alpha > 0)) & ((dist_mat_beta < treshold_beta) & (dist_mat_beta > 0)))

        d = {'r1': matrix_cutoff[0], 'r2': matrix_cutoff[1]}
        df_net = pd.DataFrame(data=d)
        df_net.name = clonotypes.name
        return df_net


