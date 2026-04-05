import numpy as np
from scipy.stats import wasserstein_distance
from immune_nets.analysis.statDistanceTypes import PairwiseStatDistance

from immune_nets.entities import GraphStats

class WassersteinStatDistance(PairwiseStatDistance):
    def __init__(self, distributionName):
        self.distributionName = distributionName
    
    def stat_dist(self, a:GraphStats, b:GraphStats) -> float:
        a_distribution: dict = getattr(a, self.distributionName)
        b_distribution: dict = getattr(b, self.distributionName)
        return wasserstein_distance(u_values=list(a_distribution.keys()),
                                    u_weights=list(a_distribution.values()),
                                    v_values=list(b_distribution.keys()),
                                    v_weights=list(b_distribution.values()),
                                   )
        


