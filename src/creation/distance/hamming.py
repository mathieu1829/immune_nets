from scipy.spatial.distance import hamming


class hammingDistance:
    def tcr_dist(self,x,ref):
        return -hamming(list(x), list(ref))
    def __str__(self):
        return "hamming"
