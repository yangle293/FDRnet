#!/usr/bin/python

from __future__ import print_function
import networkx as nx
import numpy as np
import math
import sys, argparse

# Parse arguments.
def get_parser():
    description = 'Approximately calculate the local graph with size K around a given seed.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d','--seed',type=int,required=True,help='Seed index')
    parser.add_argument('-a','--alpha',type=float,required=False,default=0.85,help='Random walk parameter')
    parser.add_argument('-s','--size',type=int,required=False,default=5,help='Local graph size')
    return parser



# Compute approximate Personalized PageRank scores given a seed.
# Adaptively change the parameter \epsilon to only calculate the first K (local graph size) largest scores

def approximatePPR_size_fast(graph, alpha, seed, size):
    isinstance(graph, nx.Graph)
    maxiter = 1e6
    epsilon = 1e-5
    not_converged = False
    d = np.array([1.0*val for (node, val) in graph.degree()])
    r = np.array([0.0]*len(d))
    p = np.array([0.0]*len(d))
    seed_index = list(nx.nodes(graph)).index(seed)
    r[seed_index] = 1
    resvec = r/d
    W = nx.adjacency_matrix(graph)
    
    I = np.where(resvec > epsilon)[0]
    # if the number of non-zero nodes is not large enough
    while len(np.where(p > 0.0)[0]) < size:
        # adaptatively change the \epsilon
        if len(np.where(p > 0.0)[0]) != 0:
            step = size/float(len(np.where(p > 0.0)[0]))
        else:
            step = 2
        epsilon = epsilon/step
        iter_n = 0
        I = np.where(resvec > epsilon)[0]
        # push the scores to neighbors
        while (I.size != 0) and (iter_n < maxiter):
            k = len(I)
            iter_n = iter_n + k
            for c in I:
                n_iter = math.ceil(math.log(r[c] / (epsilon * d[c])) / math.log(2./(1-alpha))+np.finfo(float).eps)
                p, r = push_n(pp=p, res=r, uu=c, alpha_ppr=alpha, W = W,d=d, n=n_iter)
                resvec = r/d
                I = np.where(resvec > epsilon)[0]
                ind = np.argsort(resvec[I])
                I = I[ind]

    p_vector = dict(zip(nx.nodes(graph),p))
#    if iter_n > maxiter:
#        not_converged = True
    return p_vector, not_converged, r

# push operation in PageRank-Nibble algorithm
def push_n(pp, res, uu, alpha_ppr, W,d, n):
    mult = (1 - ((1. - alpha_ppr)/2)**n) / (1 - (1 - alpha_ppr)/2)
    pp[uu] = pp[uu] + alpha_ppr * res[uu] * mult
    res = res + (1-alpha_ppr)*res[uu]/(2*d[uu])*mult*np.squeeze(W[:,uu].toarray())
    res[uu] = (1 - alpha_ppr)**n * res[uu] / (2**n)
    return pp, res

# testing function
def test(args):
    G = nx.Graph()
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 4)
    G.add_edge(4, 5)
    G.add_edge(5, 6)
    G.add_edge(6, 7)
    G.add_edge(6, 7)
    G.add_edge(7, 8)
    G.add_edge(8, 9)
    p, flag, r = approximatePPR_size_fast(G,alpha=args.alpha,seed=args.seed,size=args.size)
    print(p)

if __name__ == "__main__":
    test(get_parser().parse_args(sys.argv[1:]))
