import numpy as np
import pandas as pd

def morisita_horn_similarity(x_stats,y_stats):
    #x/y_stats are placeholders for the time being
    #it will have to be specified whether simpson index should be from tcrb or tcra or both
    #the same with "num of unique"
    #the same with "unique_distribution"
    return float( 2 * (np.array([ x_stats.unique_distribution[i] * y_stats.unique_distribution[i] for i in x_stats.num_of_unique])).sum() ) / float( (x_stats.simpson_index + y_stats.simpson_index) * x_stats.num_of_unique * y_stats.num_of_unique )

def alpha_similarity(stat_list,q):
    return float( 1 / ( ( (np.array([ [ for ] for repertoire in stat_list]).sum() )**(1/q) ) )
