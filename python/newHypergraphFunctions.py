#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Some additional functions for computing on hypergraphs


"""

import numpy as np
import random
from hypergraph import *
import scipy.stats as stats


def read_edge_data(filelocation):
    '''
    Read the data from the supplied .txt, returning it as a list of tuples. 
    '''
    with open(filelocation, 'r') as f:
        F = f.read()
    L = F.split('\n')
    L = [tuple([int(v) for v in g.split(' ') if v is not '']) for g in L]
    return L

def representing_graph(C):
    '''
    Compute the representing graph of a hypergraph, i.e., the graph with all hyperedges replaced by cliques of 2-edges.
    '''
    
    G = nx.Graph() # empty graph
    G.add_nodes_from(C.nodes) # nodes are the same as in the hypergraph
    
    for hyperedge in C.C: # go over all edges
         G.add_edges_from(itertools.combinations(hyperedge, 2))
    return(G)


def randomHypergraph(C):
    '''
    Random hypergraph that keeps only the number of hyperedges of each cardinaility
    The Erdos-Renyi-type hypergraph
    '''

    # 1) Create a random edge list
    hyperEdgesList = list()
    indexThisEdge=0
    for hyperedge in C.C: # go over all edges
        edgeCardinality = len(C.C[indexThisEdge]) # number of nodes in this hyperedge
        randomNodes = random.choices(C.nodes, k=edgeCardinality) # choose random without replacement
        hyperEdgesList.append(randomNodes) # add this hyperedge to list
        indexThisEdge = indexThisEdge+1
    # 2) Construct the hypergraph
    R = hypergraph(hyperEdgesList)
#    if C.node_labels is not None:
#        R = hypergraph(hyperEdgesList,node_labels = C.node_labels)
#    else:
#        R = hypergraph(hyperEdgesList)
            
    return(R)


def hierarchyHypergraph(hypergraph):
    # Compute the hierarchy coefficient of a hypergraph
    
    # compute degree and local clsutering coefficient
    degrees = hypergraph.D
    localClusteringThisHypergraph = localClusteringHypergraph(hypergraph)
    
    # log of degree and clustering (ignore those with degree zero because their clustering is not defined)
    
    #nonzeroDegree  = ~ np.isnan(localClusteringThisHypergraph)
    print(degrees)
    print(localClusteringThisHypergraph)
    
    nonzeroDegree  = np.where(degrees >0)
    
    logHypergraphDegree = np.log10(degrees[nonzeroDegree])
    logHypergraphClustering= np.log10(localClusteringThisHypergraph[nonzeroDegree])
    
    
    correlationResult = stats.pearsonr(logHypergraphDegree,logHypergraphClustering)
    

    return(-correlationResult[0]) 


def returnNeighbours(hypergraph,node):
    # return the neighbours of a single node in a hypergraph
    # returns a TUPLE of nodes
    
    neighbours = tuple() # empty set
    # go over each hyperedge the node is incident to
    for hyperedge in hypergraph.get_edges(node):
        # add them (elements already in there wont be added again)
       neighbours = neighbours + hyperedge
    return(tuple(set(neighbours))) # this removes all duplicate entries
    


def returnAllNeighbours(hypergraph,nodeSet):
    # retuns the set of neighbours of a set of nodes
    # returns a TUPLE of nodes
    
    neighbours = tuple() # empty set
    
    # go over each node and get all neighbours
    for node in list(nodeSet):
        # add them (elements already in there wont be added again)
       neighbours = neighbours + returnNeighbours(hypergraph,node)
    return(tuple(set(neighbours))) # this removes all duplicate entries
    

# This is the straightforward implementation. It would be faster (but use more memory) to save the extra overlap for all pairs of edges first
# or alternatively to save the extra overlap for pairs of edges that have been already computed and then just use this.
def localClusteringHypergraph(hypergraphIn):
    # local clustering for a hypergraph
    nodes = hypergraphIn.nodes # vector with nodes
    nNodes = len(nodes) # number of nodes
    
    localClustering = np.zeros(nNodes) # output for the local clustering
    
    # go over each node
    for i in range(nNodes):
        edgesThisNode = hypergraphIn.get_edges(i) # edges incident to this node
        degreeThisNode = len(edgesThisNode)
        if (degreeThisNode==0):
            localClustering[i] = float('nan')
        elif (degreeThisNode==1):
            localClustering[i] = 0
        else:
            totalEdgeOverlap = 0
            # go over all pairs of edges pairwise
            
            for edgeNumber1 in range(degreeThisNode):
                e1 = edgesThisNode[edgeNumber1]
                for edgeNumber2 in range(edgeNumber1,degreeThisNode):
                    e2 = edgesThisNode[edgeNumber2]
                    # set differences for the hyperedges
                    D1 = set(e1) - set(e2)
                    D2 = set(e2) - set(e1)
                    # if edges are the same by definition the extra overlap is zero
                    if len(D1.union(D2))==0:
                        extraOverlap=0
                    else:
                    # otherwise we have to look at their neighbours    
                        # the neighbours of D1 and D2, respectively.
                        neighD1 = set(returnAllNeighbours(hypergraphIn,tuple(D1)))
                        neighD2 = set(returnAllNeighbours(hypergraphIn,tuple(D2)))
                        # compute the extra overlap [len() is used for cardinality of edges]
                        extraOverlap = (len( neighD1.intersection(D2)  )  + len(neighD2.intersection(D1) )) /len(D1.union(D2)) # add it up
                    # add it up
                    totalEdgeOverlap = totalEdgeOverlap   +  extraOverlap
            
            # include normalisation by degree k*(k-1)/2
            localClustering[i]=2*totalEdgeOverlap/(degreeThisNode*(degreeThisNode-1))
        
    return(localClustering)

