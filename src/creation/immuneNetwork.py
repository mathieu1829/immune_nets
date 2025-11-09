import pandas as pd
import uuid

class ImmuneNetwork:
    def __init__(self,
                 graph,
                 proportions,
                 method,
                 distanceFun,
                 sampleSize,
                 sampleId,
                 threshold=0.0,
                 threshold_alpha=0.0,
                 threshold_beta=0.0,
                 name=None,
                 network_id=uuid.uuid4()
                 ):
        assert type(graph) == type(pd.DataFrame()), f"Net, should be a DataFrame, and not {type(graph)}" 
        # assert len(graph.columns) == 2, f"Net must have two columns"
        self.graph = graph 
        self.method = method
        self.sampleId = sampleId
        self.network_id = network_id
        self.distanceFun = distanceFun 
        self.threshold = threshold
        self.sampleSize = sampleSize
        self.proportions = proportions
        self.threshold_alpha = threshold_alpha
        self.threshold_beta = threshold_beta
        self.name = name
    

