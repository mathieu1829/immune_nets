#return information on graphlets within the the network for comparizon
import igraph as ig

from immune_nets.creation.algorithms.common_methods import *
from immune_nets.creation.algorithms.simpleDistance import *
from immune_nets.creation.enums.matrices import *
from immune_nets.creation.enums.utils import * 
from .immuneNetwork import ImmuneNetwork
from immune_nets.creation.utils.pathManager import pathManager


path = pathManager().testDataPath / "test_clonotypes.csv"



class GraphStats:
    def __init__(self,
                 immuneNet: ImmuneNetwork
                ):
        self.numEdges = immuneNet.graph.shape[0]
        self.nonIsolatedVerticeArray = np.unique(immuneNet.graph.to_numpy().flatten())
        self.verticeNum = self.nonIsolatedVerticeArray.shape[0]
        isolatedVerticeArray = np.setdiff1d(np.arange(immuneNet.sampleSize),self.nonIsolatedVerticeArray)
        self.isolatedVerticeNum = isolatedVerticeArray.shape[0]
        self.isolatedVerticeRatio = self.isolatedVerticeNum / self.verticeNum if self.verticeNum != 0 else -1

        self.graph = ig.Graph(immuneNet.graph.to_numpy())
        self.graph.add_vertices(immuneNet.sampleSize - self.graph.vcount())

        self.edgeDensity = float(self.graph.ecount()) / float( 0.5 * self.verticeNum * (self.verticeNum-1) ) if self.numEdges > 0 else 0.0
        self.density = self.graph.density()
        self.eccentricity = min(self.graph.eccentricity(self.nonIsolatedVerticeArray)) if self.nonIsolatedVerticeArray.size != 0 else 0 

        self.degreeDistribution = { leftBound:float(num/immuneNet.sampleSize) for leftBound, _, num in self.graph.degree_distribution().bins()}
        self.meanDegree = self.graph.degree_distribution().mean

        self.components = self.graph.components()
        self.componentSizeArray = np.array([ len(i) for i in self.components])
        self.componentCount = self.componentSizeArray.shape[0]
        self.giantComponentSize = self.components.giant().vcount()
        self.componentSizeDistribution = { component_size:(float((self.componentSizeArray == component_size).sum())/float(self.componentCount)) for component_size in np.unique(self.componentSizeArray)}
        self.proportionCountDistribution = { proportion:(immuneNet.proportions == proportion).sum() for proportion in np.unique(immuneNet.proportions)}
        componentMembers = { componentSize:[] for componentSize in np.unique(self.componentSizeArray)}
        for component in self.components:
            componentMembers[len(component)].extend(component)
        self.componentProportionDistribution = {key:immuneNet.proportions[component].sum() for key, component in componentMembers.items()}


    def toStatVector(self):
        return [ 
                float(self.isolatedVerticeRatio),
                float(self.density),
                float(self.eccentricity),
                float(self.giantComponentSize),
                float(self.meanDegree),
                float(self.componentCount)
               ]

    def toList(self):
        return [ 
                self.isolatedVerticeRatio,
                self.edgeDensity,
                self.density,
                self.eccentricity,
                self.giantComponentSize,
                self.degreeDistribution,
                self.componentCount,
                self.componentSizeDistribution,
                self.proportionCountDistribution,
                self.componentProportionDistribution,
               ]


    @staticmethod
    def vectorStatNames():
        return [
                "isolatedVerticeRatio",
                "density",
                "eccentricity",
                "giantComponent",
                "meanDegree",
                "componentCount",
               ]
    @staticmethod
    def vectorStatNamesPolish():
        return [
                "Proporcja izolowanych wierzchołków",
                "Gęstość",
                "Ekscentryczność",
                "Wielki komponent",
                "Średni stopień",
                "Liczka komponentów",
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
                "proportionCountDistribution",
                "componentProportionDistribution",
               ]
