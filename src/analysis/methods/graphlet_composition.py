#return information on graphlets within the the network for comparizon
import igraph as ig
import numpy as np
import pandas as pd

def graphletComposition(immuneNet):
    edges = immuneNet.network.shape[0]
    vertices = np.unique(immuneNet.network.to_numpy().flatten())
    vertice_num = vertices.shape[0]
    isolated_vertices = [ i for i in np.arange(immuneNet.sampleSize) if not i in vertices ]

    #transform
    # active = isolated_vertices
    # outOfBound = vertices[vertices > vertice_num]
    # minGraph = immuneNet.network.to_numpy()
    # for i in outOfBound:
    #     np.place(minGraph, minGraph == outOfBound, active.pop(0))



    # graph = ig.Graph(minGraph)
    graph = ig.Graph(immuneNet.network.to_numpy())
    graph.add_vertices(immuneNet.sampleSize - graph.vcount())

    edge_density = float(graph.ecount()) / float( 0.5 * vertice_num * (vertice_num-1) )
    percolation_threshold = immuneNet.threshold
    density = graph.density()
    eccentrity = graph.eccentrity()
    eigenvector_centrality = graph.eigenvector_centrality()
    harmonic_centrality = graph.harmonic_centrality()
    giant_component = graph.components().giant().vcount()
    betweenness = graph.betweenness()
    diameter = graph.diameter()
    closeness = graph.closeness()
    assortativity = graph.assortativity()
    assortativity_degree = graph.assortativity_degree()
    
    

    paths = np.array(graph.distances(vertices))
    mean_shortest_path = paths[paths != float('inf')].mean()

    pagerank_distribution = np.array(graph.pagerank()) 


    degree_distribution = np.array(graph.degree_distribution()) 

    components = graph.components()
    componentList = np.array([ len(i) for i in components])
    component_count = componentList.shape[0]
    component_size_distribution = { component_size:(float((componentList == component_size).sum())/float(component_count)) for component_size in np.unique(componentList)}
    return [ vertice_num,isolated_vertices,edge_density,percolation_threshold,density,eccentrity,eigenvector_centrality,harmonic_centrality,giant_component,betweenness,diameter,closeness, assortativity, assortativity_degree, mean_shortest_path, pagerank_distribution, degree_distribution, component_count,component_size_distribution]
