setExpressionData<-function(){
library(BioNet)
library(DLBCL)
library(genefilter)
library(impute)
library(igraph)
library(R.matlab)
library(locfdr)
# Overall Object: generate data
# 1. adjacent matrix of network
# 2. weighted adjacent matrix 
# 3. p-values

# load data
data(exprLym)
data(interactome)

# lym-induced network
network <- subNetwork(featureNames(exprLym), interactome) #lym-induced subnetwork
network <- largestComp(network)
network_igraph = igraph.from.graphNEL(network)
edge_list <- as_edgelist(network_igraph)
expressions <- impute.knn(exprs(exprLym))$data
weight = numeric(length=ecount(network_igraph))

# compute pearson coefficient as edge weight
for (i in c(1:ecount(network_igraph)))
{
  gene1 = edge_list[i,1]
  gene2 = edge_list[i,2]
  expr1 = expressions[gene1,]
  expr2 = expressions[gene2,]
  corr = cor(expr1,expr2,method="pearson") #absolute pearson correlation
  weight[i] = abs(corr)
}

# set weight
weighted_network = set_edge_attr(network_igraph,"weight",index=E(network_igraph),weight)
adj = as.matrix(as_adjacency_matrix(weighted_network,names=FALSE))
weight_adj = as.matrix(as_adjacency_matrix(weighted_network,attr="weight",names=FALSE))
node = V(weighted_network)$name

## t-test
t.test <- rowttests(expressions, fac = exprLym$Subgroup)
pval <- t.test$p.value
expr_name = rownames(t.test)
## locfdr and q-value
zz = qnorm(pval)
w = locfdr(zz)
lfdr = w$fdr
q = p.adjust(pval,method="fdr")

## save data
save.image(file='../othermethods/BioNet_expr/Expr_data.RData')
writeMat('../othermethods/jActiveModule_expr/Expr_data.mat',adj=adj,node=node,weight_adj=weight_adj,pval=pval,lfdr=lfdr,expr_name=expr_name,expressions=expressions,subgroup=as.numeric(exprLym$Subgroup),q=q)
writeMat('../othermethods/clusterEx_expr/Expr_data.mat',adj=adj,node=node,weight_adj=weight_adj,pval=pval,lfdr=lfdr,expr_name=expr_name,expressions=expressions,subgroup=as.numeric(exprLym$Subgroup),q=q)
writeMat('../othermethods/RegMOD_expr/Expr_data.mat',adj=adj,node=node,weight_adj=weight_adj,pval=pval,lfdr=lfdr,expr_name=expr_name,expressions=expressions,subgroup=as.numeric(exprLym$Subgroup),q=q)
writeMat('../othermethods/WMAXC_expr/Expr_data.mat',adj=adj,node=node,weight_adj=weight_adj,pval=pval,lfdr=lfdr,expr_name=expr_name,expressions=expressions,subgroup=as.numeric(exprLym$Subgroup),q=q)
writeMat('../HPC_expr/Expr_data.mat',adj=adj,node=node,weight_adj=weight_adj,pval=pval,lfdr=lfdr,expr_name=expr_name,expressions=expressions,subgroup=as.numeric(exprLym$Subgroup),q=q)
writeMat('../HotNet2_expr/Expr_data.mat',adj=adj,node=node,weight_adj=weight_adj,pval=pval,lfdr=lfdr,expr_name=expr_name,expressions=expressions,subgroup=as.numeric(exprLym$Subgroup),q=q)
}




# plot_degree_distribution = function(graph) {
#   # calculate degree
#   d = degree(graph, mode = "all")
#   dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
#   degree = 1:max(d)
#   probability = dd[-1]
#   # delete blank values
#   nonzero.position = which(probability != 0)
#   probability = probability[nonzero.position]
#   degree = degree[nonzero.position]
#   # plot
#   plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
#        col = 1, main = "Degree Distribution")
# }
# 
# plot_degree_distribution(network_igraph)
# 
# fit_power_law = function(graph) {
#   # calculate degree
#   d = degree(graph, mode = "all")
#   dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
#   degree = 1:max(d)
#   probability = dd[-1]
#   # delete blank values
#   nonzero.position = which(probability != 0)
#   probability = probability[nonzero.position]
#   degree = degree[nonzero.position]
#   reg = lm(log(probability) ~ log(degree))
#   cozf = coef(reg)
#   power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
#   alpha = -cozf[[2]]
#   R.square = summary(reg)$r.squared
#   print(paste("Alpha =", round(alpha, 3)))
#   print(paste("R square =", round(R.square, 3)))
#   # plot
#   plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
#        col = 1, main = "Degree Distribution")
#   curve(power.law.fit, col = "red", add = T, n = length(d))
# }
# 
# fit_power_law(network_igraph)
# 
# expressions <- impute.knn(exprs(exprLym))$data
# pval <- t.test$p.value
# names(pval) = names(t.test)
# fb <- fitBumModel(pval, plot = FALSE)
# scores <- scoreNodes(network = network, fb = fb, fdr = 0.001)
# network <- rmSelfLoops(network)
# module <- runFastHeinz(network, scores)
# 
# plotModule(module,scores=scores)
