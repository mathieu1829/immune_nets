import igraph as ig
import matplotlib.pyplot as plt
from pathlib import Path
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.hamming import hammingDistance
from src.creation.algorithms.simple_beta_distance import *
from src.creation.algorithms.simple_distance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.io_strategies.test_csv_strategy import *
from src.creation.immuneRepertoire import immuneRepertoire


# Example: Create a sample graph
# path = "/home/myc0plasmus/Downloads/T_PLL_sorted_5pv2_nextgem_vdj_t_clonotypes.csv"
# path = "/home/myc0plasmus/Downloads/T_PLL_sorted_5pv2_HT_nextgem_vdj_t_clonotypes.csv"
path = "/home/myc0plasmus/Downloads/10k_BMMNC_5pv2_nextgem_intron_10k_BMMNC_5pv2_nextgem_intron_vdj_t_clonotypes.csv"
# path = "/home/myc0plasmus/Documents/python/immune_nets/tests/test_data/bigTest.csv"
repertoire = test_csv_strategy().input(path)
immuneNet = simple_beta_distance(repertoire=repertoire, distance=sequenceAligner("BLOSUM62"), threshold=0.2)

print("Network:")
print(immuneNet.network.to_numpy())
edges = immuneNet.network.shape[0]
print(f"Num of edges: {edges}")
vertices = np.unique(immuneNet.network.to_numpy().flatten())
print(f"vertices: {vertices}")
prop = repertoire.clones["proportion"].to_numpy()[vertices]
print(f"proportions: {prop}")
print(f"proportions: {prop.sum()/repertoire.clones['proportion'].to_numpy().sum()}")
vertice_num = vertices.shape[0]
print(f"Num vertices: {vertice_num}")
isolated_vertices = [ i for i in np.arange(immuneNet.sampleSize) if not i in vertices]
print("Isolated vertices:")
print(isolated_vertices)
print(f"all vertices: {immuneNet.sampleSize}")


# graph = ig.Graph(minGraph)
graph = ig.Graph(immuneNet.network.to_numpy())
graph.add_vertices(immuneNet.sampleSize - graph.vcount())


# Optional: Set vertex labels or other attributes
graph.vs["label"] = [str(i) for i in range(graph.vcount())]
graph.vs["color"] = "skyblue"
bins = [0.1 * i for i in range(1, 10)]
graph.vs["size"] = (np.digitize(repertoire.clones['proportion'].to_numpy(), bins) + 3)**1.5
graph.es["width"] = 1

# Plotting the graph
layout = graph.layout("fr") # Options: "fr" (Fruchterman-Reingold), "kk", "circle", etc.

fig, ax = plt.subplots(figsize=(16, 16))
ig.plot(
    graph,
    target=ax,
    layout=layout,
    vertex_label=None,
    vertex_color=graph.vs["color"],
    vertex_size=graph.vs["size"],
    edge_width=graph.es["width"],
)
plt.show()
