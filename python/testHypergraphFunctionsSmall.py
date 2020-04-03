#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Testing the new hypergraph functions on the email network


"""


from hypergraph import * # library by P Chodrow
import numpy as np
import seaborn as sns
import pickle

from newHypergraphFunctions import * # new functions


# load the hypergraph




L = read_edge_data('./../data/smallExampleHypergraph.txt') # read the data
smalllHypergraph = hypergraph(L)     # construct a hypergraph from the data


# construct an ER-like random hypergraph
ERlikeRandomHypergraphSmall = randomHypergraph(smalllHypergraph)

# compute degree 
smalllHypergraphDegree = smalllHypergraph.D
# compute local clustering (this can be slow for large graphs)
smallHypergraphLocalClustering= localClusteringHypergraph(smalllHypergraph)






