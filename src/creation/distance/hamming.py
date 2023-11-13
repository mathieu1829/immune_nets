from scipy.spatial.distance import hamming

def hammingDistance(x,ref):
    return hamming(list(x), list(ref))
