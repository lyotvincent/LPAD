#!/usr/bin/env python
# -*- coding: utf-8 -*-

from utils import reader, normalizer
import numpy as np


class PyRWR:
    normalized = False

    def __init__(self):
        pass

    def read_graph(self, input_path, graph_type):
        '''
        Read a graph from the edge list at input_path

        inputs
            input_path : str
                path for the graph data
            graph_type : str
                type of graph {'directed', 'undirected', 'bipartite'}
        '''

        # self.A, self.base = reader.read_graph(input_path, graph_type)
        # print(self.A.shape)
        self.A = input_path
        self.base = 0
        self.m, self.n = self.A.shape
        self.node_ids = np.arange(0, self.n) + self.base
        self.normalize()

    def normalize(self):
        '''
        Perform row-normalization of the adjacency matrix
        '''
        if self.normalized is False:
            nA = normalizer.row_normalize(self.A)
            self.nAT = nA.T
            self.normalized = True
