import igraph as ig
import matplotlib.pyplot as plt
from pathlib import Path
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.hamming import hammingDistance
from src.creation.algorithms.simpleBetaDistance import *
from src.creation.algorithms.simpleDistance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.entities import ImmuneRepertoire
from src.factories import ImmuneRepertoireFactory

def graphVisualization(immuneNet, ax):

    ### DEBUG
    # print("Network:")
    # print(immuneNet.network.to_numpy())
    # edges = immuneNet.network.shape[0]
    # print(f"Num of edges: {edges}")
    # vertices = np.unique(immuneNet.network.to_numpy().flatten())
    # print(f"vertices: {vertices}")
    # prop = repertoire.clones["proportion"].to_numpy()[vertices]
    # print(f"proportions: {prop}")
    # print(f"proportions: {prop.sum()/repertoire.clones['proportion'].to_numpy().sum()}")
    # vertice_num = vertices.shape[0]
    # print(f"Num vertices: {vertice_num}")
    # isolated_vertices = [ i for i in np.arange(immuneNet.sampleSize) if not i in vertices]
    # print("Isolated vertices:")
    # print(isolated_vertices)
    # print(f"all vertices: {immuneNet.sampleSize}")

    # Creating igraph object from immuneNet
    graph = ig.Graph(immuneNet.graph.to_numpy())
    graph.add_vertices(immuneNet.sampleSize - graph.vcount())


    # Setting vertex attributes
    graph.vs["label"] = [str(i) for i in range(graph.vcount())]
    graph.vs["color"] = "skyblue"
    bins = [0.1 * i for i in range(1, 10)]
    graph.vs["size"] = (np.digitize(immuneNet.proportions, bins) + 3)**1.5
    graph.es["width"] = 1

    # Plotting the graph
    layout = graph.layout("fr") # Options: "fr" (Fruchterman-Reingold), "kk", "circle", etc.
    ig.plot(
        graph,
        target=ax,
        layout=layout,
        vertex_label=None,
        vertex_color=graph.vs["color"],
        vertex_size=graph.vs["size"],
        edge_width=graph.es["width"],
    )








               








