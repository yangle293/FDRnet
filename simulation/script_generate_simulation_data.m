%%%%% Script for generating simulation data
% Author: Le Yang
% This script is used to generated synthetic data for
% subnetwork detection in the paper:
%   FDRnet: A Novel Method to Identify Signicantly 
%   Perturbed Subnetworks in Cancer


clear all;close all;clc
%% 1. load network
[index_start,index_end] = import_edge_list('irefindex9_edge_list');
[index,gene] = import_gene_index('irefindex9_index_gene');
G = graph(gene(index_start),gene(index_end));
%% 2. Decide used complexes
% load protein complexes database CORUM
load coreComplexes.mat 
complexes = coreComplexes.subunitsGenename;
complexes = cellstr(complexes);
cc = cellfun(@(x) split(x,';'), complexes,'un',0);
cc = cellfun(@(x) upper(x),cc,'un',0);
cc_name = cellstr(coreComplexes.ComplexName);
cc_ID = coreComplexes.ComplexID;
% map to PPI and calculate sizes
for i = 1:length(cc)
    tmp = cc{i};
    tmp_network = intersect(G.Nodes.Name,tmp); %map to PPI
    if isempty(tmp_network)
        cc{i} = [];
        continue;
    end
    h = G.subgraph(tmp_network);
    aa = conncomp(h); tt = tabulate(aa);[~,indx] = max(tt(:,2)); % extract the largest component
    gene = h.Nodes.Name(find(aa==tt(indx,1)));
    cc{i} = gene;
end
cc_length = cellfun(@(x) length(x),cc,'un',0);
cc_length = cell2mat(cc_length);  % size of complexes
% load Cancer related complexes
% Paper: Detection of dysregulated protein-association networks 
%           by high-throughput proteomics predicts cancer vulnerabilities
cancer_complex = readtable('Cancer_complex_from_CORUM.xlsx');
cancer_complex_name = cancer_complex.ComplexDescription;
% human + cc_length>=10 + cc_length <= 50 + cancer related
[~,~,IB] = intersect(cancer_complex_name,coreComplexes.ComplexName);
human_index = find(strcmp(cellstr(coreComplexes.Organism),'Human'));
indx = intersect(find(cc_length>=10),human_index);
ind = intersect(indx,IB);
ind = intersect(ind,find(cc_length<=50));
% maxima non-overlapped set
ind_final = ind([1,2,3,4,5,10,12,13,17,18,19,22,25,28,29,31]); 
cc_use = cc(ind_final);
cc_name_use = cc_name(ind_final);
cc_genes = vertcat(cc_use{:});
cc_length = cc_length(ind_final);
cc_ID = cc_ID(ind_final);
usedComplex = table(cc_ID,cc_name_use,cc_length,'VariableNames',{'ID','Name','Size'});
writetable(usedComplex,'usedComplex.xls')
save('usedComplexes.mat','cc_use')
%% 3. Assign p-values using mixture model: beta+uniform
[~,sig_index,~] = intersect(G.Nodes.Name,cc_genes);
nonsig_index = setdiff(1:length(G.Nodes.Name),sig_index);

n_permute = 10;
beta_a = [0.11,0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01];
p_all = cell(length(beta_a),n_permute);
q_all = cell(length(beta_a),n_permute);
z_all = cell(length(beta_a),n_permute);
for i = 1:length(beta_a)
    for j = 11:n_permute
        rng('shuffle')
        score = zeros(length(G.Nodes.Name),1);
        score(nonsig_index) = rand([length(nonsig_index),1]);
        score(sig_index) = random('beta',beta_a(i),1,length(sig_index),1);
        q = mafdr(score,'bhfdr',true);
        q = -log10(q);
        zz = -norminv(score);
        p_all{i,j} = score;
        q_all{i,j} = q;
        z_all{i,j} = zz;
%         table_hotnet = table(G.Nodes.Name,q);
         table_FDRnet = table(G.Nodes.Name,zz);
%         writetable(table_hotnet,['hotdata/simu_q_beta_',num2str(beta_a(i)),'_',num2str(j),'.txt'],'WriteVariableNames',0,'Delimiter','\t');
         writetable(table_FDRnet,['fdrdata/simu_z_beta_',num2str(beta_a(i)),'_',num2str(j),'.txt'],'WriteVariableNames',0,'Delimiter','\t');       
    end
end
%save('Data.mat','p_all','q_all','z_all')
save('Data_0.01.mat','p_all','q_all','z_all')

