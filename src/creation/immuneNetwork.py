import pandas as pd

class immuneNetwork:
    def __init__(self,network,method, sampleIDs, distanceFun, threshold, sampleSize):
        assert type(network) == type(pd.DataFrame()), f"Net, should be a DataFrame, and not {type(network)}" 
        assert len(network.columns) == 2, f"Net must have two columns"
        self.algorithmParams = {}
        self.graph = network 
        self.method = method
        self.sampleIDs = sampleIDs
        self.distanceFun = distanceFun 
        self.algorithmParams["threshold"] = threshold
        self.sample_size = sampleSize
    

