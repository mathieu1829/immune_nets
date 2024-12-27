from tests.utils.cdr3SeqGeneratorUtil import *
import numpy as np

class cdr3SetGenerator:
    def generateRandomSubset(self, n, min_length=11, max_length=16):
        repertoireData = [ cdr3SeqGeneratorUtil().generate_random_cdr3(min_length=min_length,max_length=max_length) for _ in range(n)]
        return np.array(repertoireData)

    #To do - create a more deterministic version
    def generateRelatedSubset(self, n, distance, distanceFun):
        centralSequence = cdr3SeqGeneratorUtil().generate_random_cdr3()
        subset = np.array([cdr3SeqGeneratorUtil().mutateCDR3(centralSequence, numMutations=random.randint(1, 3)) for _ in range(n)])
        # assert isinstance(centralSequence, str) and len(centralSequence) > 0, "Central sequence is invalid."
        # assert all(isinstance(seq, str) and len(seq) > 0 for seq in subset), "One or more sequences in the subset are invalid."
        distances = np.array([distanceFun(centralSequence, seq) for seq in subset])
        iterCount = 0
        while((distances > distance).sum() != 0):
            correctedSubset = np.array([cdr3SeqGeneratorUtil().mutateCDR3(centralSequence) for _ in range((distances > distance).sum()) ])
            subset[distances > distance] = correctedSubset
            distances = np.array([distanceFun(centralSequence, seq) for seq in subset])
            iterCount += 1
            if(iterCount > 100):
                raise ValueError("Correcting this subset took more than one hundred iteration, are you sure it can be corrected with provided parameters?")
        
        return subset

