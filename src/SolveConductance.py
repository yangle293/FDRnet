from __future__ import print_function
import cplex
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import product
from cplex.exceptions import CplexError, CplexSolverError
import time


def solveConductance(G_undirected, bound_fdr, terminal, full_Graph, time_limit, relative_gap):

    ############## Build ILP problem #############
    ## Variables include:
    # Z: Objective function
    # Zx: helper variables --> nNodes
    # xx: helper variables --> nNodes**2
    # binary indicators: x --> nNodes,
    # flow on edges: y_edge --> nEdges,
    # residual flow: z_0 --> 1
    # source to node edges: y_source to terminal --> 1
    # number of nodes and edges in network
    
    ## Basic numbers
    G = nx.DiGraph(G_undirected)
    nNodes = G.number_of_nodes()
    nEdges = G.number_of_edges()
    node_list = list(nx.nodes(G))
    edge_list = nx.edges(G)
    
    #full_laplacian = nx.laplacian_matrix(full_Graph) # Laplacian matrix for full graph D - A
    #full_node_list = list(nx.nodes(full_Graph)) # full node list
    #index = [full_node_list.index(x) for x in node_list] # index of the small graph
    #laplacian = full_laplacian[index,:] # extract Laplacian matrix of small graph
    #laplacian = laplacian[:,index]
    
    locfdr = nx.get_node_attributes(G, 'scores') # get local fdr dict
    full_deg_list = dict(full_Graph.degree()) # degree dict for full graph
    deg_list = [(node,full_deg_list[node]) for node in node_list] # degree dict for small graph
    adj = nx.adjacency_matrix(G);
    
    u = 1.001 # upper bound of Z
    l = -0.001 # lower bound of Z

    ######## Variables Start
    # Name
    variable_x = [x for x in nx.nodes(G)] # x_i
    variable_z = "obj" # Z
    variable_zx = ["".join(("z", x)) for x in nx.nodes(G)] # zx_i
    variable_xx = ["".join((x, "*", y)) for (x,y) in nx.edges(G)] # xx_ij #edges
    variable_y = ["".join((x, "_", y)) for x, y in nx.edges(G)] #y_ij #egdes
    variable_source = "source" # source flow
    variable_s2x = ["".join((variable_source, "_", terminal))] # source to terminal
    my_colnames = variable_x + [variable_z] + variable_zx + variable_xx +  variable_y + [variable_source] + variable_s2x
    variable_dict = dict(zip(my_colnames,range(len(my_colnames)))) # variable map, speed up problem building of cplex

    # Type: x_i binary
    my_ctype = "I" * nNodes + "C" * (nEdges + nNodes + nEdges + 3)
    # bound: x {0,1}, obj 0-1, zx: 0-1, xx: 0-1, y: 0-nodes, source: 0-nodes, source2seed: 0-nodes
    my_ub = [1.0] * nNodes + [1.0] * (1 + nNodes + nEdges) + [float(nNodes)] * (nEdges + 2)
    my_lb = [0.0] * nNodes + [0.0] * (1 + nNodes + nEdges + nEdges + 2)
    ######## Variables End
    
    ######## Objective
    my_obj = [0.0] * nNodes + [1.0] + [0.0] * (nEdges + nNodes + nEdges + 2) # objective: min Z
    
    ######## Constraints Start
    my_row = []
    #constrain_start = time.time()
    
    #### Terminal constrain  x_terminal = 1
    my_row.append([[variable_dict[terminal]], [1.0]])
    my_rhs_terminal = [1.0]
    my_sense_terminal = "E"
    my_rownames_terminal = ["terminal"]
    
    #### Budget constraint fdr: sum(lfdr-B) <= 0
    my_row.append([ [variable_dict[x] for x in variable_x], [ float(locfdr[x] - bound_fdr) for x in variable_x if x in locfdr.keys()] ] )
    my_sense_budget_fdr = "L"
    my_rhs_budget_fdr = [0.0]
    my_rownames_budget_fdr = ["Budget_fdr"]
    
    #### Connectivity constraints
    
    # residual flow + flow injected into network = total flow  z0 + sum(y_zv) = N
    # my_row.append([["z0"]+my_colnames_s2x, [1.0]*(nNodes+1)])
    my_row.append([[variable_dict[variable_source], variable_dict[variable_s2x[0]]], [1.0] * (1 + 1)])
    my_sense_total = "E"
    my_rhs_total = [float(nNodes)]
    my_rownames_total = ["total"]
    
    # positive flow
    for i in range(nEdges):
        my_row.append([[variable_dict[variable_y[i]], variable_dict[variable_y[i].split("_")[-1]]], [1.0, -float(nNodes)]])
    my_rhs_positive = [0.0] * nEdges
    my_sense_positive = "L" * nEdges
    my_rownames_positive = ["".join((variable_y[i], "flow")) for i in range(nEdges)]

    # consuming flow
    my_rhs_consuming = []
    my_sense_consuming = "E" * nNodes
    my_rownames_consuming = []
    for node in nx.nodes(G):
        edge_incoming = [variable_dict["".join((x,"_",y))] for (x, y) in G.in_edges(node)]
        edge_ongoing = [variable_dict["".join((x,"_",y))] for (x, y) in G.out_edges(node)]
        if node == terminal:
            edge_incoming.append(variable_dict[variable_s2x[0]])
        edge_ongoing.append(variable_dict[node])
        my_row.append([edge_incoming + edge_ongoing, [1.0] * len(edge_incoming) + [-1.0] * (len(edge_ongoing))])
        my_rhs_consuming = my_rhs_consuming + [0.0]
        my_rownames_consuming = my_rownames_consuming + ["".join(("consuming_", node))]

    # injected = consumed
    my_rhs_equal = [0.0]
    my_sense_equal = "E"
    my_rownames_equal = ["Equal"]
    # my_row.append([my_colnames_x + my_colnames_s2x, [1.0]*nNodes+[-1.0]*nNodes])
    tmp = [variable_dict[x] for x in variable_x]
    tmp.append(variable_dict[variable_s2x[0]])
    my_row.append([tmp, [1.0] * nNodes + [-1.0]])

#print("finish connectivity constraints...")
    #### Helper constraints
    # Main: -sum(zx_i * D_ii) + sum((D_ij-A_ij)*x_ij) <= 0
    # New Main: -sum(zx_i * D_ii) + sum(D_ii*x_ii) - sum(A_ij*x_ij) <= 0 and x_ii == x_i
    my_rhs_obj = [0.0]
    my_sense_obj = "L"
    my_rownames_obj = ["obj_main"]
    var_zx = [variable_dict[x] for x in variable_zx]
    var_x = [variable_dict[x] for x in variable_x]
    var_xx = [variable_dict[x] for x in variable_xx]
    tmpvar = []
    tmpvar.extend(var_zx + var_x + var_xx)

    value_zx = [-1.0*val for (node,val) in deg_list]
    value_x = [1.0*val for (node,val) in deg_list]
    value_xx = [-1.0*adj[node_list.index(x.split('*')[0]),node_list.index(x.split('*')[1])] for x in variable_xx]
    tmpvalue = []
    tmpvalue.extend(value_zx + value_x + value_xx)

#    tmp1 = [laplacian[i,j] for (i,j) in product(range(len(node_list)),range(len(node_list)))]
#    tmp2 = [-1.0*val for (node,val) in deg_list]
#    tmp_var = [variable_dict[x] for x in variable_xx]
#    tmp_var.extend([variable_dict[x] for x in variable_zx])
#    tmp1.extend(tmp2)
    my_row.append([tmpvar,tmpvalue])


    # Constraints for zx:
    # 1. zx_i >= z - u(1-x_i) ===>   z +   x_i - zx_i <= u
    # 2. zx_i >= l*x_i        ===>       l*x_i - zx_i <= 0
    # 3. zx_i <= z - l(1-x_i) ===>  -z - l*x_i + zx_i <= -l
    # 4. zx_i <= u*x_i        ===>     - u*x_i + zx_i <= 0
    my_rhs_zx = [u] * nNodes + [0.0] * nNodes + [-1.0*l] * nNodes + [0.0] * nNodes
    my_sense_zx = "L" * 4 * nNodes
    my_rownames_zx = ["".join(("obj_helper_zx1_",node_list[i])) for i in range(len(node_list))] + ["".join(("obj_helper_zx2_",node_list[i])) for i in range(len(node_list))] + ["".join(("obj_helper_zx3_",node_list[i])) for i in range(len(node_list))] + ["".join(("obj_helper_zx4_",node_list[i])) for i in range(len(node_list))]

    constrain_zx_1 = [[[variable_dict[variable_z],variable_dict[variable_x[i]],variable_dict[variable_zx[i]]],[1.0,u,-1.0]] for i in range(len(node_list))]
    constrain_zx_2 = [[[variable_dict[variable_x[i]],variable_dict[variable_zx[i]]],[l,-1.0]] for i in range(len(node_list))]
    constrain_zx_3 = [[[variable_dict[variable_z],variable_dict[variable_x[i]],variable_dict[variable_zx[i]]],[-1.0,-l,1.0]] for i in range(len(node_list))]
    constrain_zx_4 = [[[variable_dict[variable_x[i]],variable_dict[variable_zx[i]]],[-u,1.0]] for i in range(len(node_list))]
    my_row.extend(constrain_zx_1 + constrain_zx_2 + constrain_zx_3 + constrain_zx_4)
#print("finishing zx constraints...")


    # Constraints for xx
    # 1. xx_ij <= x_i                ===>  xx_ij - x_i        <= 0
    # 2 xx_ij <= x_j               ===>  xx_ij - x_j        <= 0
    # 3 xx_ij >= x_i + x_j - 1     ===> -xx_ij + x_i + x_j  <= 1

    my_rhs_xx = [0.0] * (2 * nEdges) + [1.0] * nEdges
    my_sense_xx = "L" * (3 * nEdges)
    my_rownames_xx = ["".join(("obj_helper_xx1_",x.split('*')[0],x.split('*')[1])) for x in variable_xx] + \
                    ["".join(("obj_helper_xx2_",x.split('*')[0],x.split('*')[1])) for x in variable_xx] + \
                    ["".join(("obj_helper_xx3_",x.split('*')[0],x.split('*')[1])) for x in variable_xx]
#my_rownames_xx = ["".join(("obj_helper_xx1_",node_list[i],node_list[j])) for (i,j) in product(range(len(node_list)),range(len(node_list))) if i!=j] + ["".join(("obj_helper_xx2_",node_list[i],node_list[j])) for (i,j) in product(range(len(node_list)),range(len(node_list))) if i!=j] + ["".join(("obj_helper_xx2_",node_list[i],node_list[i])) for i in range(len(node_list))] + ["".join(("obj_helper_xx3_",node_list[i],node_list[j])) for (i,j) in product(range(len(node_list)),range(len(node_list))) if i!=j] + ["".join(("obj_helper_xx3_",node_list[i],node_list[i])) for i in range(len(node_list))]

    constrain_1 = [[[variable_dict[x],variable_dict[x.split('*')[0]]],[1.0,-1.0]] for x in variable_xx]
    constrain_2 = [[[variable_dict[x],variable_dict[x.split('*')[1]]],[1.0,-1.0]] for x in variable_xx]
    constrain_3 = [[[variable_dict[x],variable_dict[x.split('*')[0]],variable_dict[x.split('*')[1]]],[-1.0,1.0,1.0]] for x in variable_xx]
#    constrain_3_2 = [[[variable_dict["".join((node_list[i], "*", node_list[i]))],variable_dict[variable_x[i]]],[-1.0,2.0]] for i in range(len(node_list))]
    my_row.extend(constrain_1 + constrain_2 + constrain_3)
#    print("finishing xx constraints...")
#constrain_end = time.time()
#    print("constrain time: ", constrain_end - constrain_start)

    #### Merge Constraints
    my_sense = my_sense_terminal + my_sense_budget_fdr + my_sense_total + my_sense_positive + my_sense_consuming + my_sense_equal + \
                my_sense_obj + my_sense_zx + my_sense_xx
    my_rhs = my_rhs_terminal + my_rhs_budget_fdr + my_rhs_total + my_rhs_positive + my_rhs_consuming + my_rhs_equal + \
                my_rhs_obj + my_rhs_zx + my_rhs_xx
    my_rownames = my_rownames_terminal + my_rownames_budget_fdr + my_rownames_total + \
                    my_rownames_positive + my_rownames_consuming + my_rownames_equal + my_rownames_obj + my_rownames_zx + my_rownames_xx
#    print(variable_dict)
#    for (x,y,z) in zip(my_rownames,my_row,my_rhs):
#        print(x,y,z)
    ######## Constraints End
    
    ######## Problem Build Start
    try:
        build_start = time.time()
        subnetwork = cplex.Cplex()
        subnetwork.parameters.timelimit.set(time_limit)
        subnetwork.parameters.mip.tolerances.mipgap.set(relative_gap)
        subnetwork.set_log_stream(None)
        subnetwork.set_warning_stream(None)
        subnetwork.set_results_stream(None)
        #subnetwork.parameters.simplex.tolerances.markowitz.set(0.99999)
        #subnetwork.parameters.read.scale.set(1)
        subnetwork.objective.set_sense(subnetwork.objective.sense.minimize)
        subnetwork.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, types=my_ctype, names=my_colnames)
        subnetwork.linear_constraints.add(lin_expr=my_row, senses=my_sense,rhs=my_rhs, names=my_rownames)
        #build_end = time.time()
        #print("build time: ", build_end - build_start)
        #solve_start = time.time()
        subnetwork.solve()
        solve_end = time.time()
#print("Solve time: ", solve_end - solve_start)
    except CplexError as exc:
        #print(exc)
        subnetwork_result = [];
        status = -1;
        return
    ######## Problem Build End
    
    ######## Solution Start
    solution = subnetwork.solution
    #solution.get_status() returns  an integer code
#print("Solution status = ", solution.get_status(), ":", end=' ')
    # # the following line prints the corresponding string
#    print(solution.status[solution.get_status()])
    status = solution.status[solution.get_status()]
    try:
        #       print("Objective value = ", solution.get_objective_value())
        x = solution.get_values(0, subnetwork.variables.get_num() - 1)
        result = dict(zip(my_colnames, x))
    except CplexSolverError as exc:
        #print(exc)
        subnetwork_result = []
        return subnetwork_result,status
#print(result)
    subnetwork_result = [x for x in variable_x if result[x] > 0.9]
    return subnetwork_result, status
    ######## Solution End
    
#print(my_colnames)
#   print(my_obj)




if __name__ == "__main__":
    G = nx.Graph()
    G.add_nodes_from(["TP53","PTEN","PIK3CA","MAP3K1"])
    G.add_edges_from([("TP53","PTEN"),("PIK3CA","MAP3K1"),("TP53","MAP3K1"),("TP53","PIK3CA"),
                      ("TP53","KKK"),("PTEN","MAP3K1"),("PIK3CA","PTEN")])
    scores = {"TP53":0.0,"PTEN":0.0,"PIK3CA":0,"MAP3K1":1.0,"KKK":1.0}
    nx.set_node_attributes(G,scores,'scores')
    locfdr = nx.get_node_attributes(G, 'scores')
    subnetwork_result,status = solveConductance(G,0.25,"TP53",G)
    print(subnetwork_result)
