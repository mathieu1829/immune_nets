import numpy as np
from scipy import stats as scipy_stats

from .distributionIdentifier import DistributionIdentifier

class KolomogorovIndentifier(DistributionIdentifier):
    def __init__(self, usePvalue=True):
        self.usePvalueFlag = usePvalue

    def identify_distribution(self, statDict):
        values = np.array(list(statDict.keys()))
        probs  = np.array(list(statDict.values()))
        if probs.sum() != 1:
            probs = probs / probs.sum()

        results = {}
        for name, dist in self.candidates.items():
            # generate samples from discrete distribution
            sample = np.random.choice(values, size=100000, p=probs)
            # try to fit the to theoretical candidate distribution
            try:
                params = dist.fit(sample)

                theoretical = dist.pdf(values, *params)
                theoretical /= theoretical.sum()
               
                usePvalue = 1 if self.usePvalueFlag else 0

                results[name] = scipy_stats.kstest(sample, name, params)[usePvalue]
                # if that fails assign high score to guarantee choosing different distribution
            except:
                results[name] = 100000
        if self.usePvalueFlag :
            best = max(results, key=results.get) # type: ignore
        else:
            best = min(results, key=results.get) # type: ignore

        
        return best
