#return information on graphlets within the the network for comparizon
import igraph as ig
import numpy as np
import pandas as pd
import math

import uuid
import src.creation.algorithms.simpleDistance 
import src.creation.distance.alignment
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simpleDistance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.utils.pathManager import pathManager
from src.factories import ImmuneRepertoireFactory
from src.mappers import GraphStatsMapper


path = pathManager().testDataPath / "test_clonotypes.csv"



class GraphStats:
    def __init__(self,immuneNet):
        edges = immuneNet.graph.shape[0]
        vertices = np.unique(immuneNet.graph.to_numpy().flatten())
        self.verticeNum = vertices.shape[0]
        self.isolatedVertices = [ i for i in np.arange(immuneNet.sampleSize) if not i in vertices ]
        self.isolatedVerticeNum = len(self.isolatedVertices)
        self.isolateVerticeRatio = self.isolatedVerticeNum / self.verticeNum

        #transform
        # active = isolated_vertices
        # outOfBound = vertices[vertices > vertice_num]
        # minGraph = immuneNet.network.to_numpy()
        # for i in outOfBound:
        #     np.place(minGraph, minGraph == outOfBound, active.pop(0))



        # graph = ig.Graph(minGraph)
        self.graph = ig.Graph(immuneNet.graph.to_numpy())
        self.graph.add_vertices(immuneNet.sampleSize - self.graph.vcount())

        self.edgeDensity = float(self.graph.ecount()) / float( 0.5 * self.verticeNum * (self.verticeNum-1) ) if edges > 0 else 0.0
        self.density = self.graph.density()
        self.eccentrity = np.array(self.graph.eccentricity())
        self.giantComponent = self.graph.components().giant().vcount()
        # assortativity = graph.assortativity()
        # assortativity_degree = graph.assortativity_degree()
        
        

        self.paths = np.array(self.graph.distances(vertices))

        self.pagerankDistribution = self.graph.pagerank() 

        self.degreeDistribution = self.graph.degree_distribution()
        self.meanDegree = self.degreeDistribution.mean

        self.components = self.graph.components()
        self.componentList = np.array([ len(i) for i in self.components])
        self.componentCount = self.componentList.shape[0]
        self.componentSizeDistribution = { component_size:(float((self.componentList == component_size).sum())/float(self.componentCount)) for component_size in np.unique(self.componentList)}
        self.meanComponentSize = sum([ key*self.componentSizeDistribution[key] for key in self.componentSizeDistribution])


if __name__ == "__main__":
    df_net = simpleDistance(repertoire=ImmuneRepertoireFactory.fromCSV(path=path, name="test repertoire", desc=""), distance=sequenceAligner("BLOSUM62"))
    graphletList = GraphStatsMapper.toList(GraphStats(df_net))
    print(graphletList)
