import os
import sys

from pyrwr.rwr import RWR
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from colormap import Color, Colormap
import argparse


def rwr_graph(input_graph,
              graph_type='undirected',
              c=0.85,
              epsilon=1e-9,
              max_iters=100,
              shape=0
              ):
    """
    Random Walk with Restart

    Params:
    input_graph - an N by N Symmetric Matrix
    graph_type - undirected | directed | bipartite
    c - restart probability
    epsilon - a small float num, see Method
    max_iters - Stop Early
    shape - the column or row of input_graph | label size

    Returns:
        return an N by N matrix P^t, see Method
    """
    rwr = RWR()
    vectors = []
    rwr.read_graph(input_path=input_graph, graph_type=graph_type)
    for i in range(shape):
        r = rwr.compute(i, c, epsilon, max_iters)
        vectors.append(r)
    vectors = np.asarray(vectors)
    return vectors


def matrix_to_graph(matrix, file, p):
    """
    process matrix --

    Params:

    matrix - the P^t
    file - the file path
    p - penalty coefficient

    Returns:
        None
    """
    f = open(file, 'w')
    shape = matrix.shape[0]
    for i in range(shape):
        for j in range(i + 1, shape):
            if matrix[i][j] == 0:
                continue
            w = float(matrix[i][j]) / abs(i - j) ** p
            line = str(i) + ' ' + str(j) + ' ' + str(w) + '\r\n'
            f.write(line)


def label_propagation(vector_dict, edge_dict):
    """
    Label Propagation Process

    Params:
    vector_dict - all node, like {str -> int}
    edge_dict -  all edges, like  {edge1：['node1:weight1','node2:weight2']}

    Returns:
        After t Round, the weight Matrix P^t
    """
    t = 0
    while True:
        if (check(vector_dict, edge_dict) == 0):
            t = t + 1
            for node in vector_dict.keys():
                adjacency_node_list = edge_dict[node]
                vector_dict[node] = get_max_community_label(
                    vector_dict, adjacency_node_list)
        else:
            break
    return vector_dict


def tadpcc(s, e, m):
    """
    Function for compute PCC value

    Params:
        s - TAD start bin index
        e - TAD end bin index
        m - contact matrix

    Returns
        return the PCC value of a TAD
    """
    res = 0
    cnt = 0
    if s > e:
        s, e = e, s
    for si in range(s, e + 1):
        for sj in range(si + 1, e + 1):
            crf = np.corrcoef(m[si], m[sj])[0][1]
            cnt += 1
            if np.isnan(crf):
                res += 1
            else:
                res += crf
    if cnt == 0:
        # print(s,e)
        return 0
    return res / cnt


def taddiff(s0, e0, s1, e1, s2, e2, m):
    """
    Function for compute DIFF

    Params:
        s0 - the upper TAD's start bin index
        e0 - the upper TAD's end bin index
        s1 - this TAD's start bin index
        e1 - this TAD's end bin index
        s2 - the down TAD's start bin index
        e2 - the down TAD's end bin index
        m - the contact matrix
    """
    res = 0
    cnt1 = 0
    cnt2 = 0
    if s0 > e0:
        s0, e0 = e0, s0
    if s1 > e1:
        s1, e1 = e1, s1
    if s2 > e2:
        s2, e2 = e2, s2
    # Intra(TAD)
    intra = 0
    for i in range(s1, e1 + 1):
        for j in range(i, e1 + 1):
            intra += m[i][j]
            cnt1 += 1
    # Inter(TAD_{i-1} , TAD_i)
    intre = 0
    # upper
    if s0 != 0:
        for i in range(s0, e0 + 1):
            for j in range(s1, e1 + 1):
                intre += m[i][j]
                cnt2 += 1
    # down
    if s2 != 0:
        for i in range(s1, e1 + 1):
            for j in range(s2, e2 + 1):
                intre += m[i][j]
                cnt2 += 1
    intra = intra / cnt1
    intre = intre / cnt2
    return intra - intre


def get_max_community_label(vector_dict, adjacency_node_list):
    """
    Function for getting all node's labels from its adjacency node

    Params:
        vector_dict - all node, like {str -> int}
        adjacency_node_list - a dict that stored a node's adjacency_node, like {node1:[node2, node3]}

    Returns:
        t_th Rounds, each node and its labels


    """
    label_dict = {}
    for node in adjacency_node_list:
        node_id_weight = node.strip().split(":")
        node_id = node_id_weight[0]  # neighbor
        # the weight between a node and its neighbor
        node_weight = float(node_id_weight[1])
        if vector_dict[node_id] not in label_dict:
            label_dict[vector_dict[node_id]] = node_weight
        else:
            label_dict[vector_dict[node_id]] += node_weight
    sort_list = sorted(label_dict.items(), key=lambda d: d[1], reverse=True)
    return sort_list[0][0]


def check(vector_dict, edge_dict):
    """
    Determine that all nodes have the correct labels
    """
    for node in vector_dict.keys():
        adjacency_node_list = edge_dict[node]  # all nodes link to this node
        node_label = vector_dict[node]  # the label of this node
        label = get_max_community_label(vector_dict, adjacency_node_list)
        if node_label == label:
            continue
        else:
            return 0
    return 1


def loadData(filePath):
    """
    create two dict from P^t

    Params:
        filePath: where P^t is stored

    Returns:
        vector_dict - all node, like {str -> int}
        edge_dict -  all edges, like  {edge1：['node1:weight1','node2:weight2']}
    """
    f = open(filePath)
    vector_dict = {}
    edge_dict = {}
    for line in f.readlines():
        lines = line.strip().split(" ")
        for i in range(2):
            if lines[i] not in vector_dict:
                vector_dict[lines[i]] = int(lines[i])
                edge_list = []
                if len(lines) == 3:
                    edge_list.append(lines[1 - i] + ":" + lines[2])
                else:
                    edge_list.append(lines[1 - i] + ":" + "1")
                edge_dict[lines[i]] = edge_list
            else:
                edge_list = edge_dict[lines[i]]
                if len(lines) == 3:
                    edge_list.append(lines[1 - i] + ":" + lines[2])
                else:
                    edge_list.append(lines[1 - i] + ":" + "1")
                edge_dict[lines[i]] = edge_list
    return vector_dict, edge_dict


def LPA(filePath):
    """
    Label propagation
    """
    vector, edge = loadData(filePath)
    vector_dict = label_propagation(vector, edge)
    cluster_group = dict()
    for node in vector_dict.keys():
        cluster_id = vector_dict[node]
        if cluster_id not in cluster_group.keys():
            cluster_group[cluster_id] = [node]
        else:
            cluster_group[cluster_id].append(node)
    return cluster_group


def _plot_HiC(matrix_data, vmax, colors=None):
    """
    Function for plot Hi-C heat-map

    Params:
        matrix_data - contact matrix
        vmax - make all records that lager vmax to be max
        colors - color bar

    Returns:
        A Figure
    """
    if colors is None:
        colors = ['white', 'red']
    red_list = list()
    green_list = list()
    blue_list = list()
    for color in colors:
        col = Color(color).rgb
        red_list.append(col[0])
        green_list.append(col[1])
        blue_list.append(col[2])
    c = Colormap()
    d = {'blue': blue_list,
         'green': green_list,
         'red': red_list}
    mycmap = c.cmap(d)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_facecolor('w')
    ax.grid(b=None)
    sns.heatmap(
        matrix_data.T,
        vmax=vmax,
        xticklabels=100,
        yticklabels=100,
        cmap=mycmap,
        cbar=False)


def com_local_density(matrix, w=5, topk=60):
    """
    compute local density

    Params:
     matrix - P^t
     w - window size
     top-k - 0 ~ 100 ,
    """
    shape = matrix.shape[0]
    local_density = np.zeros((shape,))
    mean_Val = np.mean(matrix)
    for i in range(shape):
        for l in range(i - w, i):
            for m in range(i, i + w):
                if l < 0 or m >= shape:
                    local_density[i] += mean_Val
                else:
                    local_density[i] += matrix[l][m]
    minimize = []
    for i in range(1, shape - 1):
        if local_density[i -
                         1] > local_density[i] and local_density[i] < local_density[i +
                                                                                    1]:
            minimize.append(i)
    minimize = np.asarray(minimize)
    lg = []
    rg = []
    for i in minimize:
        j = i - 1
        while j >= 0 and local_density[j] > local_density[j + 1]:
            j -= 1
        # lg.append((local_density[j + 1] - local_density[i]) / (i - j))
        lg.append(i - j)
        j = i + 1
        while j + 1 < shape and local_density[j] < local_density[j + 1]:
            j += 1
        # rg.append((local_density[j] - local_density[i]) / (j - i))
        rg.append(j - i)
    lg = np.asarray(lg)
    rg = np.asarray(rg)
    gr = lg + rg
    res = [0]
    if topk == 0:
        for i in range(len(gr)):
            if gr[i] > 2*w:
                res.append(minimize[i])
    else:
        thr = np.percentile(gr, topk)
        for i in range(len(gr)):
            if gr[i] >= thr:
                res.append(minimize[i])
    res.append(shape)
    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Detection of Topological Associated Domains from Hi-C Data using Network Construction and Label propagation')
    parser.add_argument(
        '-f',
        type=str,
        help="the path of a intra-chromosomal Hi-C matrix seperated by Tab with N by N shaped")
    parser.add_argument(
        '-w',
        type=int,
        default=5,
        help="window size,the default is 5.")
    parser.add_argument(
        '-c',
        type=float,
        default=0.9,
        help="the restart probability, the default is 0.9")
    parser.add_argument(
        '-o',
        type=str,
        default='./out',
        help="the storage path and filename of result, the default is ./out")
    parser.add_argument(
        '-p',
        type=float,
        default=1.0,
        help="optional, the Penalty coefficient, the default is 1")
    parser.add_argument(
        '-k',
        type=float,
        default=0.6,
        help="optional, the top k, 0.1 ~ 0.9， the default value is inferred based on w automatically. Not recommended to use")
    args = parser.parse_args()
    filepath, w, c, out, p, topk = args.f, args.w, args.c, args.o, args.p, args.k

    if filepath is None:
        print("no file input")
        sys.exit()

    if topk:
      topk = (1 - topk) * 100
    else:
      topk = 0

    ###########################################################################
    # Step 1: Read contact matrix from input file
    ###########################################################################
    f = open(out, 'w')
    bak_matrix = np.loadtxt(filepath)

    ###########################################################################
    # Step 2: RWR to Create Graph
    ###########################################################################
    bak_matrix = rwr_graph(bak_matrix, c=c, shape=bak_matrix.shape[0])

    ###########################################################################
    # Step 3: compute local density and select cut-offs
    ###########################################################################
    local_density = com_local_density(bak_matrix, w=w, topk=topk)
    local_density.insert(0, 0)
    local_density.append(bak_matrix.shape[0])

    ###########################################################################
    # Step 4: label propagation in each sub-graph
    ###########################################################################
    for bins in range(0, len(local_density) - 1):
        s = local_density[bins]
        e = local_density[bins + 1]
        matrix = bak_matrix[s:e, s:e]
        # matrix = rwr_graph(matrix, c=c, shape=matrix.shape[0])
        lpa_file = './matrix_for_lpa'
        matrix_to_graph(matrix, file=lpa_file, p=(p))
        cluster_g = LPA(lpa_file)
        os.remove(lpa_file)

        #######################################################################
        # Step 4.2: filter and export TAD from communities
        #######################################################################
        for g in cluster_g:
            t = cluster_g[g]
            up = int(t[0])
            down = int(t[-1])
            if down - up < 4:
                continue
            else:
                x = up + s
                y = down + s
                f.write(str(x) + '\t' + str(y) + '\n')
    f.close()
