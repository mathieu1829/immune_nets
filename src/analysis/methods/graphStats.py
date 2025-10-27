#return information on graphlets within the the network for comparizon
import igraph as ig

from src.creation.algorithms.common_methods import *
from src.creation.algorithms.simpleDistance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.immuneNetwork import ImmuneNetwork
from src.creation.utils.pathManager import pathManager


path = pathManager().testDataPath / "test_clonotypes.csv"



class GraphStats:
    def __init__(self,
                 immuneNet: ImmuneNetwork
                ):
        edges = immuneNet.graph.shape[0]
        vertices = np.unique(immuneNet.graph.to_numpy().flatten())
        self.verticeNum = vertices.shape[0]
        isolatedVertices = [ i for i in np.arange(immuneNet.sampleSize) if not i in vertices ]
        self.isolatedVerticeNum = len(isolatedVertices)
        self.isolatedVerticeRatio = self.isolatedVerticeNum / self.verticeNum

        self.graph = ig.Graph(immuneNet.graph.to_numpy())
        self.graph.add_vertices(immuneNet.sampleSize - self.graph.vcount())

        self.edgeDensity = float(self.graph.ecount()) / float( 0.5 * self.verticeNum * (self.verticeNum-1) ) if edges > 0 else 0.0
        self.density = self.graph.density()
        self.eccentricity = np.array(self.graph.eccentricity())
        self.meanEccentricity = self.eccentricity.mean()
        self.giantComponent = self.graph.components().giant().vcount()

        self.degreeDistribution = self.graph.degree_distribution()
        self.meanDegree = self.degreeDistribution.mean

        self.components = self.graph.components()
        self.componentList = np.array([ len(i) for i in self.components])
        self.componentCount = self.componentList.shape[0]
        self.componentSizeDistribution = { component_size:(float((self.componentList == component_size).sum())/float(self.componentCount)) for component_size in np.unique(self.componentList)}
        self.meanComponentSize = sum([ key*self.componentSizeDistribution[key] for key in self.componentSizeDistribution])

    def toStatVector(self):
        return [ 
                float(self.isolatedVerticeRatio),
                float(self.edgeDensity),
                float(self.density),
                float(self.meanEccentricity),
                float(self.giantComponent),
                float(self.meanDegree),
                float(self.componentCount),
                float(self.meanComponentSize)
               ]

    def toList(self):
        return [ 
                self.isolatedVerticeRatio,
                self.edgeDensity,
                self.density,
                self.eccentricity,
                self.giantComponent,
                self.degreeDistribution,
                self.componentCount,
                self.componentSizeDistribution
               ]


    @staticmethod
    def vectorStatNames():
        return [
                "isolatedVerticeRatio",
                "edgeDensity",
                "density",
                "meanEccentricity",
                "giantComponent",
                "meanDegree",
                "componentCount",
                "meanComponentSize"
               ]
    @staticmethod
    def listStatNames():
        return [
                "isolatedVerticeRatio",
                "edgeDensity",
                "density",
                "eccentricity",
                "giantComponent",
                "degreeDistribution",
                "componentCount",
                "componentSizeDistribution"
               ]


