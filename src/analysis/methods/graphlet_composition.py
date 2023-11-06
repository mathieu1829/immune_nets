#return information on graphlets within the the network for comparizon
import igraph as ig
import numpy as np
import pandas as pd

def graphletComposition(raw_network):
    edges = raw_network.shape[0]
    vertices = np.unique(raw_network.to_numpy().flatten())
    vertice_num = vertices.shape[0]
    graph = ig.GraphBase(vertice_num, raw_network.to_numpy())
    edge_density = float(graph.ecount()) / float( 0.5 * vertice_num * (vertice_num-1) )
    #percolation_threshold = ummm threshold?
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

    pagerank_of_nonexisting_vertices = np.delete(np.array(graph.pagerank()),vertices).sum()
    pagerank_distribution = np.array([ i + (pagerank_of_nonexisting_vertices / vertices.shape[0]) for i in np.array(graph.pagerank())[vertices]])
   #pagerank_distribution = np.array(graph.pagerank()) #or maybe just like that


    degree_distribution_of_nonexisting_vertices = np.delete(np.array(graph.degree_distribution()),vertices).sum()
    degree_distribution = np.array([ i + (degree_distribution_of_nonexisting_vertices / vertices.shape[0]) for i in np.array(graph.degree_distribution())[vertices]])
    #degree_distribution = np.array(graph.degree_distribution()) #or maybe just like that

    components = graph.components()
    componentList = np.array([ components[i].vcount for i in range(components) if components[i].vcount != 1 ])
    component_count = componentList.shape[0]
    component_size_distribution = { component_size:(float(count)/float(component_count)) for component_size, count in np.unique(componentList)}
    # component_count = len(graph.components())
    # componentList = np.array([ components[i].vcount for i in range(components) ])
    # component_size_distribution = { component_size:(float(count)/float(component_count)) for component_size, count in np.unique(componentList)}
