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
        self.numEdges = immuneNet.graph.shape[0]
        self.nonIsolatedVertices = np.unique(immuneNet.graph.to_numpy().flatten())
        self.verticeNum = self.nonIsolatedVertices.shape[0]
        isolatedVertices = [ i for i in np.arange(immuneNet.sampleSize) if not i in self.nonIsolatedVertices ]
        self.isolatedVerticeNum = len(isolatedVertices)
        self.isolatedVerticeRatio = self.isolatedVerticeNum / self.verticeNum if self.verticeNum != 0 else -1

        self.graph = ig.Graph(immuneNet.graph.to_numpy())
        self.graph.add_vertices(immuneNet.sampleSize - self.graph.vcount())

        self.edgeDensity = float(self.graph.ecount()) / float( 0.5 * self.verticeNum * (self.verticeNum-1) ) if self.numEdges > 0 else 0.0
        self.density = self.graph.density()
        self.eccentricity = min(self.graph.eccentricity(self.nonIsolatedVertices)) if self.nonIsolatedVertices.size != 0 else 0 
        self.giantComponent = self.graph.components().giant().vcount()

        self.degreeDistribution = { leftBound:float(num/immuneNet.sampleSize) for leftBound, _, num in self.graph.degree_distribution().bins()}
        self.meanDegree = self.graph.degree_distribution().mean

        self.components = self.graph.components()
        self.componentList = np.array([ len(i) for i in self.components])
        self.componentCount = self.componentList.shape[0]
        self.componentSizeDistribution = { component_size:(float((self.componentList == component_size).sum())/float(self.componentCount)) for component_size in np.unique(self.componentList)}
        self.componentProportionDistribution = { idx:immuneNet.proportions[component].sum for idx,component in enumerate(self.components)}


    def toStatVector(self):
        return [ 
                float(self.isolatedVerticeRatio),
                float(self.edgeDensity),
                float(self.density),
                float(self.eccentricity),
                float(self.giantComponent),
                float(self.meanDegree),
                float(self.componentCount)
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
                self.componentSizeDistribution,
                self.componentProportionDistribution
               ]


    @staticmethod
    def vectorStatNames():
        return [
                "isolatedVerticeRatio",
                "edgeDensity",
                "density",
                "eccentricity",
                "giantComponent",
                "meanDegree",
                "componentCount",
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
                "componentSizeDistribution",
                "componentProportionDistribution"
               ]


