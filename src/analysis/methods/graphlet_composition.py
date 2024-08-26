#return information on graphlets within the the network for comparizon
import igraph as ig
import numpy as np
import pandas as pd

import uuid
import src.creation.io_strategies.df_strategy 
import src.creation.algorithms.simple_distance 
import src.creation.distance.alignment
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simple_distance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.io_strategies.test_csv_strategy import *
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.utils.pathManager import pathManager
import pickle 

path = pathManager().testDataPath / "test_clonotypes.csv"



class graphletComposition:
    def __init__(self,immuneNet):
        edges = immuneNet.network.shape[0]
        vertices = np.unique(immuneNet.network.to_numpy().flatten())
        self.vertice_num = vertices.shape[0]
        self.isolated_vertices = [ i for i in np.arange(immuneNet.sampleSize) if not i in vertices ]

        #transform
        # active = isolated_vertices
        # outOfBound = vertices[vertices > vertice_num]
        # minGraph = immuneNet.network.to_numpy()
        # for i in outOfBound:
        #     np.place(minGraph, minGraph == outOfBound, active.pop(0))



        # graph = ig.Graph(minGraph)
        self.graph = ig.Graph(immuneNet.network.to_numpy())
        self.graph.add_vertices(immuneNet.sampleSize - self.graph.vcount())

        self.edge_density = float(self.graph.ecount()) / float( 0.5 * self.vertice_num * (self.vertice_num-1) )
        self.percolation_threshold = immuneNet.threshold
        self.density = self.graph.density()
        self.eccentrity = self.graph.eccentricity()
        self.eigenvector_centrality = [round(i,6) for i in self.graph.eigenvector_centrality()]
        self.harmonic_centrality = self.graph.harmonic_centrality()
        self.giant_component = self.graph.components().giant().vcount()
        self.betweenness = self.graph.betweenness()
        self.diameter = self.graph.diameter()
        self.closeness = self.graph.closeness()
        # assortativity = graph.assortativity()
        # assortativity_degree = graph.assortativity_degree()
        
        

        self.paths = np.array(self.graph.distances(vertices))
        self.mean_shortest_path = self.paths[self.paths != float('inf')].mean()

        self.pagerank_distribution = np.array(self.graph.pagerank()) 


        self.degree_distribution = np.array(self.graph.degree_distribution()) 

        self.components = self.graph.components()
        self.componentList = np.array([ len(i) for i in self.components])
        self.component_count = self.componentList.shape[0]
        self.component_size_distribution = { component_size:(float((self.componentList == component_size).sum())/float(self.component_count)) for component_size in np.unique(self.componentList)}
    # return [ vertice_num,isolated_vertices,edge_density,percolation_threshold,density,eccentrity,eigenvector_centrality,harmonic_centrality,giant_component,betweenness,diameter,closeness, assortativity, assortativity_degree, mean_shortest_path, pagerank_distribution, degree_distribution, component_count,component_size_distribution]
    def toList(self):
        return [ self.vertice_num,self.isolated_vertices,self.edge_density,self.percolation_threshold,self.density,self.eccentrity,self.eigenvector_centrality,self.harmonic_centrality,self.giant_component,self.betweenness,self.diameter,self.closeness, self.mean_shortest_path, self.pagerank_distribution, self.degree_distribution, self.component_count,self.component_size_distribution]

if __name__ == "__main__":
    df_net = simple_distance(repertoire=test_csv_strategy().input(path), distance=sequenceAligner("BLOSUM62"))
    graphletList = graphletComposition(df_net).toList()
    print(graphletList)
