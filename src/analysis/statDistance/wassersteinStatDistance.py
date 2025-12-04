import numpy as np
from scipy.stats import wasserstein_distance

from src.analysis.methods.graphStats import GraphStats

class WassersteinStatDistance:
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
        


