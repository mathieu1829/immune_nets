import igraph as ig

import matplotlib.pyplot as plt
from pathlib import Path
from immune_nets.creation.algorithms.common_methods import *
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.creation.distance.hamming import hammingDistance
from immune_nets.creation.algorithms.simpleBetaDistance import *
from immune_nets.creation.algorithms.simpleDistance import *
from immune_nets.creation.enums.matrices import *
from immune_nets.creation.enums.utils import * 
from immune_nets.entities import ImmuneRepertoire
from immune_nets.factories import ImmuneRepertoireFactory
from matplotlib.collections import LineCollection

def graphVisualization(immuneNet, ax):

    ### debug
    # print("network:")
    # print(immunenet.network.to_numpy())
    # edges = immunenet.network.shape[0]
    # print(f"num of edges: {edges}")
    # vertices = np.unique(immunenet.network.to_numpy().flatten())
    # print(f"vertices: {vertices}")
    # prop = repertoire.clones["proportion"].to_numpy()[vertices]
    # print(f"proportions: {prop}")
    # print(f"proportions: {prop.sum()/repertoire.clones['proportion'].to_numpy().sum()}")
    # vertice_num = vertices.shape[0]
    # print(f"num vertices: {vertice_num}")
    # isolated_vertices = [ i for i in np.arange(immunenet.samplesize) if not i in vertices]
    # print("isolated vertices:")
    # print(isolated_vertices)
    # print(f"all vertices: {immunenet.samplesize}")

    # creating igraph object from immunenet
    edges = immuneNet.graph.to_numpy()  # (E, 2)

    graph = ig.Graph(
        n=immuneNet.sampleSize,
        edges=edges,
        directed=False
    )

    layout = graph.layout_fruchterman_reingold(
        niter=400,   # default ~1000 (too slow)
        # grid=True
    )  # or lgl/fr
    print("create layout")
    coords = np.asarray(layout.coords)

    bins = [0.1 * i for i in range(1, 10)]
    sizes = (np.digitize(immuneNet.proportions, bins) + 3)**1.5

    # ---- FAST EDGE DRAWING ----
    edge_coords = coords[edges]  # (E, 2, 2)
    lc = LineCollection(edge_coords, linewidths=0.3, alpha=0.3)
    ax.add_collection(lc)

    # ---- nodes ----
    ax.scatter(coords[:,0], coords[:,1], s=sizes)

    ax.axis("off")

    del graph
    del layout








               








