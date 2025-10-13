import pandas as pd
import uuid

class ImmuneNetwork:
    def __init__(self,
                 graph,
                 method,
                 distanceFun,
                 sampleSize,
                 threshold=0.0,
                 threshold_alpha=0.0,
                 threshold_beta=0.0,
                 name=None,
                 sampleId=uuid.uuid4
                 ):
        assert type(graph) == type(pd.DataFrame()), f"Net, should be a DataFrame, and not {type(graph)}" 
        assert len(graph.columns) == 2, f"Net must have two columns"
        self.graph = graph 
        self.method = method
        self.sampleId = sampleId
        self.distanceFun = distanceFun 
        self.threshold = threshold
        self.sampleSize = sampleSize
        self.threshold_alpha = threshold_alpha
        self.threshold_beta = threshold_beta
        self.name = name
    

