#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Testing the new hypergraph functions


"""


from hypergraph import * # library by P Chodrow
import numpy as np
import seaborn as sns
import pickle

from newHypergraphFunctions import * # new functions


# load the hypergraph
filehandler = open("./../data/hypergraphProteins.out", 'rb') 
proteinHypergraph = pickle.load(filehandler)

# construct an ER-like random hypergraph
ERlikeRandomHypergraph = randomHypergraph(proteinHypergraph)

# compute degree 
proteinHypergraphDegree = proteinHypergraph.D
# compute local clustering (this can be slow for large graphs)
proteinHypergraphLocalClustering= localClusteringHypergraph(proteinHypergraph)


