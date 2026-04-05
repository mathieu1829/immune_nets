import numpy as np

from .distributionIdentifier import DistributionIdentifier

class KlDivergenceIndentifier(DistributionIdentifier):
    def kl_div(self, P, Q, eps=1e-12):
        P = P + eps
        Q = Q + eps
        P /= P.sum()
        Q /= Q.sum()
        return np.sum(P * np.log(P / Q))

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

                results[name] = self.kl_div(probs, theoretical)
                # if that fails assign high score to guarantee choosing different distribution
            except:
                results[name] = 10
        best = min(results, key=results.get) # type: ignore
        
        return best
