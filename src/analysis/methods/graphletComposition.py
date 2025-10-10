#return information on graphlets within the the network for comparizon
import igraph as ig
import numpy as np
import pandas as pd
import math

import uuid
import src.creation.io_strategies.df_strategy 
import src.creation.algorithms.simpleDistance 
import src.creation.distance.alignment
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simpleDistance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.utils.pathManager import pathManager


path = pathManager().testDataPath / "test_clonotypes.csv"



class graphletComposition:
    def __init__(self,immuneNet):
        edges = immuneNet.graph.shape[0]
        vertices = np.unique(immuneNet.graph.to_numpy().flatten())
        self.vertice_num = vertices.shape[0]
        self.isolated_vertices = [ i for i in np.arange(immuneNet.sample_size) if not i in vertices ]
        self.isolated_vertices_num = len(self.isolated_vertices)

        #transform
        # active = isolated_vertices
        # outOfBound = vertices[vertices > vertice_num]
        # minGraph = immuneNet.network.to_numpy()
        # for i in outOfBound:
        #     np.place(minGraph, minGraph == outOfBound, active.pop(0))



        # graph = ig.Graph(minGraph)
        print(immuneNet.graph)
        self.graph = ig.Graph(immuneNet.graph.to_numpy())
        self.graph.add_vertices(immuneNet.sample_size - self.graph.vcount())

        self.edge_density = float(self.graph.ecount()) / float( 0.5 * self.vertice_num * (self.vertice_num-1) ) if edges > 0 else 0.0
        self.percolation_threshold = immuneNet.algorithmParams["threshold"]
        self.density = self.graph.density()
        self.eccentrity = np.array(self.graph.eccentricity())
        self.eigenvector_centrality = np.array([round(i,6) for i in self.graph.eigenvector_centrality()])
        self.harmonic_centrality = np.array(self.graph.harmonic_centrality())
        self.giant_component = self.graph.components().giant().vcount()
        self.betweenness = np.array(self.graph.betweenness())
        self.diameter = self.graph.diameter()
        self.closeness = self.graph.closeness()
        self.mean_closeness = np.array([x if not math.isnan(x) else 0 for x in self.closeness]).mean()
        # assortativity = graph.assortativity()
        # assortativity_degree = graph.assortativity_degree()
        
        

        self.paths = np.array(self.graph.distances(vertices))
        self.mean_shortest_path = self.paths[self.paths != float('inf')].mean() if edges > 0 else -1.0

        self.pagerank_distribution = self.graph.pagerank() 
        self.expected_pagerank = np.array(self.pagerank_distribution).mean()


        self.degree_distribution = self.graph.degree_distribution()
        self.expected_degree = self.degree_distribution.mean

        self.components = self.graph.components()
        self.componentList = np.array([ len(i) for i in self.components])
        self.component_count = self.componentList.shape[0]
        self.component_size_distribution = { component_size:(float((self.componentList == component_size).sum())/float(self.component_count)) for component_size in np.unique(self.componentList)}
        self.expected_component_size = sum([ key*self.component_size_distribution[key] for key in self.component_size_distribution])
    # return [ vertice_num,isolated_vertices,edge_density,percolation_threshold,density,eccentrity,eigenvector_centrality,harmonic_centrality,giant_component,betweenness,diameter,closeness, assortativity, assortativity_degree, mean_shortest_path, pagerank_distribution, degree_distribution, component_count,component_size_distribution]
    def toList(self):
        return [ 
                float(self.vertice_num),
                float(self.isolated_vertices_num),
                float(self.edge_density),
                float(self.percolation_threshold),
                float(self.density),
                float(self.eccentrity.mean()),
                float(self.eigenvector_centrality.mean()),
                float(self.harmonic_centrality.mean()),
                float(self.giant_component),
                float(self.betweenness.mean()),
                float(self.diameter),
                float(self.mean_closeness),
                float(self.mean_shortest_path),
                float(self.expected_pagerank),
                float(self.expected_degree),
                float(self.component_count),
                float(self.expected_component_size)
                ]

# if __name__ == "__main__":
#     df_net = simpleDistance(repertoire=test_csv_strategy().input(path), distance=sequenceAligner("BLOSUM62"))
#     graphletList = graphletComposition(df_net).toList()
#     print(graphletList)
