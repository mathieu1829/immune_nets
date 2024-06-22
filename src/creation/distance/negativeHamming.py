from scipy.spatial.distance import hamming

class negativeHammingDistance:
    def tcr_dist(self,x,ref):
        maxLen = max(len(x), len(ref))
        x = x.ljust(maxLen)
        ref = ref.ljust(maxLen)
        return hamming(list(x), list(ref))

    def __str__(self):
        return "negative_hamming"

