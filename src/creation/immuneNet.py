import pandas as pd

class immuneNet:
    def __init__(self,net, sampleIDs, distanceFun, threshold):
        assert type(net) == type(pd.DataFrame()), f"Net, should be a DataFrame, and not {type(net)}" 
        assert len(net.columns) == 2, f"Net must have two columns"
        self.net = net 
        self.sampleIDs = sampleIDs
        self.distanceFun = distanceFun 
        self.threshold = threshold
    

