import pandas as pd

from src.creation.immuneNetwork import ImmuneNetwork
from src.models import Network

class NetworkMapper:
    @staticmethod
    def toImmuneNetwork(network: Network):
        new_graph = pd.DataFrame([{
            "r1": n.r1,
            "r2": n.r2
        } for n in network.network_edges])

        return ImmuneNetwork(graph=new_graph,
                             method=network.algorithm,
                             sampleId=network.repertoire_id,
                             distanceFun=network.distance_function,
                             threshold=eval(network.network_algorithm_parameters)["threshold"],
                             sampleSize=len(network.source_repertoire.clonotypes),
                             name=network.name
                            )

    @staticmethod
    def fromImmuneNetwork( network: ImmuneNetwork):
        parameters =  {"threshold":network.threshold}
        parameters = str(parameters)
        new_network = Network(repertoire_id=network.sampleId,
                          algorithm=network.method,
                          distance_function=network.distanceFun,
                          network_algorithm_parameters=parameters,
                          sample_size=network.sampleSize,
                          name = network.name
                          )
        new_network.setGraph(network.graph)
        return new_network
