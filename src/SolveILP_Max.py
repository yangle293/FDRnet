from __future__ import print_function
import cplex
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from cplex.exceptions import CplexError, CplexSolverError


def solveNetworkProblem(G_undirected, bound_fdr, terminal):

    ############## Build ILP problem #############
    # Variables include:
    # binary indicators: x --> nNodes,
    # flow on edges: y_edge --> nEdges,
    # residual flow: z_0 --> 1
    # source to node edges: y_source to terminal --> 1
    # number of nodes and edges in network
    G = nx.DiGraph(G_undirected)
    nNodes = G.number_of_nodes()
    nEdges = G.number_of_edges()

    # transformed local fdr value for each node
    locfdr = nx.get_node_attributes(G, 'scores')
    coe = nx.get_node_attributes(G, 'ppr')
    #tmp = nx.get_node_attributes(G, 'ppr')
    #coe = dict((x, np.log10(y)-np.min(np.log10(tmp.values()))) for (x, y) in tmp.items())

    # variables name
    source = "source"
    my_colnames_x = [x for x in nx.nodes(G)]
    my_colnames_y = ["".join((x, "_", y)) for x, y in
                     nx.edges(G)]
    my_colnames_s2x = ["".join((source, "_", terminal))]
    my_colnames = my_colnames_x + my_colnames_y + my_colnames_s2x + [source]

    # variables type
    my_ctype = "I" * nNodes + "C" * (nEdges + 2)
    # objective function
    my_obj = [1.0] * nNodes + [0.0] * (nEdges + 2)

    # lower and upper bound
    my_ub = [1.0] * nNodes + [cplex.infinity] * (nEdges + 1) + [float(nNodes)]
    my_lb = [0.0] * (nNodes + nEdges + 2)

    # Constraints:
    my_row = []

    # terminal version
    my_row.append([[terminal], [1.0]])
    my_rhs_terminal = [1.0]
    my_sense_terminal = "E"
    my_rownames_terminal = ["terminal"]

    # budget constraint fdr
    my_row.append([my_colnames_x, [locfdr[x]-bound_fdr for x in my_colnames_x if x in locfdr.keys()]])
    my_sense_budget_fdr = "L"
    my_rhs_budget_fdr = [0.0]
    my_rownames_budget_fdr = ["Budget_fdr"]

    # residual flow + flow injected into network = total flow  z0 + sum(y_zv) = N
    # my_row.append([["z0"]+my_colnames_s2x, [1.0]*(nNodes+1)])
    my_row.append([[source, my_colnames_s2x[0]], [1.0] * (1 + 1)])

    my_sense_total = "E"
    my_rhs_total = [float(nNodes)]
    my_rownames_total = ["total"]

    # positive flow
    for i in range(nEdges):
        my_row.append([[my_colnames_y[i], my_colnames_y[i].split("_")[-1]], [1.0, -float(nNodes)]])

    my_rhs_positive = [0.0] * nEdges
    my_sense_positive = "L" * nEdges
    my_rownames_positive = ["".join((my_colnames_y[i], "flow")) for i in range(nEdges)]



    # consuming flow
    my_rhs_consuming = []
    my_sense_consuming = "E" * nNodes
    my_rownames_consuming = []
    for node in nx.nodes(G):
        edge_incoming = ["".join((x,"_",y)) for (x, y) in G.in_edges(node)]
        edge_ongoing = ["".join((x,"_",y)) for (x, y) in G.out_edges(node)]
        if node == terminal:
            edge_incoming.append(my_colnames_s2x[0])
        edge_ongoing.append(node)
        my_row.append([edge_incoming + edge_ongoing, [1.0] * len(edge_incoming) + [-1.0] * (len(edge_ongoing))])
        my_rhs_consuming = my_rhs_consuming + [0.0]
        my_rownames_consuming = my_rownames_consuming + ["".join(("consuming_", node))]

    # injected = consumed
    my_rhs_equal = [0.0]
    my_sense_equal = "E"
    my_rownames_equal = ["Equal"]
    # my_row.append([my_colnames_x + my_colnames_s2x, [1.0]*nNodes+[-1.0]*nNodes])
    tmp = list(my_colnames_x)
    tmp.append(my_colnames_s2x[0])
    my_row.append([tmp, [1.0] * nNodes + [-1.0]])
    # merge all constraints
    my_sense = my_sense_terminal + my_sense_budget_fdr + my_sense_total + my_sense_positive + my_sense_consuming + my_sense_equal
    my_rhs = my_rhs_terminal + my_rhs_budget_fdr + my_rhs_total + my_rhs_positive + my_rhs_consuming + my_rhs_equal
    my_rownames = my_rownames_terminal + my_rownames_budget_fdr + my_rownames_total + my_rownames_positive + my_rownames_consuming + my_rownames_equal
    ########### Start to build problem in CPLEX ###############
    try:
        subnetwork = cplex.Cplex()
    # parameters
        subnetwork.parameters.timelimit.set(60.0)
        #subnetwork.parameters.mip.tolerances.absmipgap.set(1e-9)
        subnetwork.set_log_stream(None)
        subnetwork.set_warning_stream(None)
        subnetwork.set_results_stream(None)
        subnetwork.objective.set_sense(subnetwork.objective.sense.maximize)
        subnetwork.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, types=my_ctype,
                       names=my_colnames)
        subnetwork.linear_constraints.add(lin_expr=my_row, senses=my_sense,
                                rhs=my_rhs, names=my_rownames)
        subnetwork.solve()
    except CplexError as exc:
        print(exc)
        return

    ###########  Get solution ##################
    solution = subnetwork.solution
    # solution.get_status() returns  an integer code
    # print("Solution status = ", solution.get_status(), ":", end=' ')
    # # the following line prints the corresponding string
    # print(solution.status[solution.get_status()])
    # print("Objective value = ", solution.get_objective_value())
    status = solution.status[solution.get_status()]
    try:
        x = solution.get_values(0, subnetwork.variables.get_num() - 1)
        result = dict(zip(my_colnames, x))
        subnetwork_result = [x for x in my_colnames_x if result[x] > 0.9]
        return subnetwork_result, status
    except CplexSolverError as exc:
        print(exc)
        subnetwork_result = []
        return subnetwork_result, status
