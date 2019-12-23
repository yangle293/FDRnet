from __future__ import print_function
import networkx as nx
import numpy as np
import sys
import os
import operator
import argparse
import matplotlib.pyplot as plt

from FDRnet_main import load_network_from_file, largest_component


# Parse arguments.
def get_parser():
    description = 'Plot FDRnet subnetworks results.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-igi','--input_gene_index',type=str,required=True,help='Input gene index file name')
    parser.add_argument('-iel','--input_edge_list',type=str,required=True,help='Input edge list file name')
    parser.add_argument('-igl','--input_gene_lfdr',type=str,required=True,help='Input gene local FDRs')
    parser.add_argument('-ofn','--output_file_name',type=str,required=True,help='File name of output')
    return parser


def plot_subnetwork(G,result):
    for r in result.keys():
        sub = G.subgraph(result[r])
        nodes = sub.nodes()
        colors = [sub.node[n]['scores'] for n in nodes]
        plt.figure(figsize=(8,5))
        pos = nx.spring_layout(sub)
        ec = nx.draw_networkx_edges(sub, pos, alpha=0.5)
        pp = nx.draw_networkx_nodes(sub, pos, alpha=0.7,node_list = nodes,node_color=colors,cmap=plt.cm.jet)
        nn = nx.draw_networkx_labels(sub,pos)
        plt.colorbar(pp)
        plt.axis('off')
        figname = "seed_"+r+".png"
        plt.savefig(figname,format="PNG")
#plt.show()
def load_subnetwork(filename):
    with open(filename,'r') as f:
        arrs = [l.rstrip().split(',') for l in f]
        arrs.pop(0)
        result = dict((arr[0],arr[3].rstrip().split()) for arr in arrs)
        return result


# Run script
def run(args):
    index_file = args.input_gene_index; edge_file = args.input_edge_list; score_file = args.input_gene_lfdr; result_file = args.output_file_name
    G_raw = load_network_from_file(index_file=index_file, edge_file=edge_file, score_file=score_file)
    G = largest_component(G_raw)
    result = load_subnetwork(result_file)
    plot_subnetwork(G,result)


if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
