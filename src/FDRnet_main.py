from __future__ import print_function
import networkx as nx
import numpy as np
import sys
import os
import operator
import heapq
import time
import argparse
import matplotlib.pyplot as plt

from APPR_fast import approximatePPR_size_fast
from SolveConductance import solveConductance
from SolveILP_Max import solveNetworkProblem

# Parse arguments.
def get_parser():
    description = 'FDRnet: identifying significant mutated subnetworks in a network.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-igi','--input_gene_index',type=str,required=True,help='Input gene index file name')
    parser.add_argument('-iel','--input_edge_list',type=str,required=True,help='Input edge list file name')
    parser.add_argument('-igl','--input_gene_lfdr',type=str,required=True,help='Input gene local FDRs')
    parser.add_argument('-ofn','--output_file_name',type=str,required=True,help='File name of output')
    parser.add_argument('-se','--seed',type=str,required=True,default='all',help='Seed gene name')
    parser.add_argument('-sm','--sort_method',type=str,required=False,default='neighbor',help='Seed sort method')
    parser.add_argument('-bd','--bound',type=float,required=False,default='0.1',help='FDR bound')
    parser.add_argument('-al','--alpha',type=float,required=False,default=0.85,help='Random walk parameter')
    parser.add_argument('-sz','--size',type=int,required=False,default=400,help='Local graph size')
    parser.add_argument('-tl','--time_limit',type=int,required=False,default=100,help='Time limit for each seed')
    parser.add_argument('-rg','--relative_gap',type=float,required=False,default=0.01,help='Relative gap in MILP')
    return parser


#####################################################################################################################################################
#
# FDRnet main functions
#
#####################################################################################################################################################


def FDRnet(attGraph, bound, alpha, size, time_limit, relative_gap, sort_method, seed):
    ## Seed selection and sort
    if seed == 'all':
        scores = nx.get_node_attributes(attGraph, 'scores') # local FDRs dictionary
        values = [(x,scores[x]) for x in attGraph.nodes()] # associate local FDRs with all genes
        seed_set = [x for (x,y) in values if y <= bound] # seeds (local FDRs less than bound)
        seed_subnetwork = attGraph.subgraph(seed_set) # seeds subnetwork
        seed_data = [(x,scores[x],seed_subnetwork.degree(x)) for x in seed_set] # (seed,lfdr,#seeds-in-neighbor)
        # sort the seeds with different criteria
        if sort_method == 'lfdr':
            seed_sort_lfdr = sorted(seed_data, key=operator.itemgetter(1))
            seed_list = [x for (x,y,z) in seed_sort_lfdr]
        elif sort_method == 'neighbor':
            seed_sort_neighbor = sorted(seed_data, key=operator.itemgetter(2),reverse=True)
            seed_list = [x for (x,y,z) in seed_sort_neighbor]
    else:
        seed_list = [seed]
    ## Start searching subnetworks
    total_start = time.time()
    selected_seed = set() # maintain a set of selected genes
    result = list()

    for seed in seed_list:
        print("Seed ", seed, str(seed_list.index(seed)+1),"/",str(len(seed_list)))
        if seed in selected_seed:
            continue
        seed_start = time.time()
        result_seed, status = denseSubgraph(attGraph, bound, alpha, size, seed, time_limit, relative_gap)
        seed_end = time.time()
        result_sub = attGraph.subgraph(result_seed)
        if result_sub.degree(seed) == 1:
            continue
        selected_seed = selected_seed.union(result_seed)
        result.append((seed, seed_end-seed_start, status, result_seed))
    total_end = time.time()
    print("Total time: ", total_end - total_start)
    return result


def denseSubgraph(attGraph, bound, alpha, size, seed, time_limit, relative_gap):
    ## TO DO: check parameters
    #
    #print("start to approximate PPR...")
    #ppr_start = time.time()
    ppr, converge, res = approximatePPR_size_fast(graph=attGraph, alpha=alpha, seed=seed, size=size)
    #ppr_end = time.time()
    #print("finished PPR computation...", ppr_end-ppr_start)
    normalized_ppr = dict((x, y/float(attGraph.degree[x])) for (x, y) in ppr.items())
    nx.set_node_attributes(attGraph, normalized_ppr, 'ppr')
    nodes = heapq.nlargest(int(size), normalized_ppr, key=normalized_ppr.get)
    G = attGraph.subgraph(nodes)
    # presolve: large problem may not get solution with size > 1 within 60s
    #print("pre-solve the problem...")
    #pre_start = time.time()
    subnetwork, status = solveNetworkProblem(G, bound, seed)
    
    #pre_end = time.time()
    #print("pre_solve time: ",pre_end-pre_start)
    if len(subnetwork) == 1 and status == 'MIP_optimal':
        #print("pre-solve find single gene as maximum solution...")
        return subnetwork, "".join(("presolve_",status))
    else:
        #print("pre-solve can not decide the result...")
        subnetwork, status = solveConductance(G, bound, seed, attGraph, time_limit, relative_gap)
        return subnetwork, status


#####################################################################################################################################################
#
# Load network data
#
#####################################################################################################################################################

def load_network_from_file(index_file, edge_file, score_file):
    # Load gene-index map
    with open(index_file) as infile:
        arrs = [l.rstrip().split() for l in infile]
        indexToGene = dict((int(arr[0]), arr[1]) for arr in arrs)

    G = nx.Graph()
    G.add_nodes_from(indexToGene.values())  # in case any nodes have degree zero

    # Load graph
    with open(edge_file) as infile:
        edges = [map(int, l.rstrip().split()[:2]) for l in infile]
    G.add_edges_from([(indexToGene[u], indexToGene[v]) for u, v in edges])

    # Load p-values or scores
    with open(score_file) as infile:
        arrs = [l.rstrip().split() for l in infile]
        geneToScores = dict((arr[0], float(arr[1])) for arr in arrs)
    # TO DO: compute local fdr from p-value
    # merge scores and set score of no-score gene as 1
    for gene in indexToGene.values():
        if gene not in geneToScores.keys():
            geneToScores[gene] = 1.0
    nx.set_node_attributes(G, geneToScores, 'scores')
    return G


# Remove self-loops, multi-edges, and restrict to the largest component
def largest_component(G):
    selfLoops = [(u, v) for u, v in G.edges() if u == v]
    G.remove_edges_from(selfLoops)
    return G.subgraph(sorted(nx.connected_components(G), key=lambda cc: len(cc), reverse=True)[0])




# Run script
def run(args):
    index_file = args.input_gene_index; edge_file = args.input_edge_list; score_file = args.input_gene_lfdr; result_file = args.output_file_name
    bound = args.bound; alpha = args.alpha; size = args.size; time_limit = args.time_limit;
    relative_gap = args.relative_gap; sort_method=args.sort_method;seed = args.seed
    
    G_raw = load_network_from_file(index_file=index_file, edge_file=edge_file, score_file=score_file)
    G = largest_component(G_raw)

    result = FDRnet(G, bound=bound, alpha=alpha, size=size,time_limit = time_limit, relative_gap = relative_gap, sort_method = sort_method, seed=seed)
#    this_folder = os.path.dirname(os.path.abspath(__file__))
#    my_file = os.path.join(this_folder, result_file)
    with open(result_file, 'wb') as outfile:
        outfile.write("".join(("Seed Gene", ",", "Running Time",",","Optimization Status",",","Subnetwork", "\n")))
        for (seed,time,status, r) in result:
            outfile.write("".join((seed, ",", str(time),",",status,","," ".join(r), "\n")))


if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))

