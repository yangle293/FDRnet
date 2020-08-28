%%% Goal: 
%%1. Use expression data to filter CNV data
%%2. Use filtered CNV data to compute p-values
%%3. Merge p-values to generate input data for all the methods
clear all;close all;clc

%%
%%%%% Parameters
nan_threshold = 0.2;

%%
%%%% Load the Expression data and Copy Number data
%load cnv
CNV_raw = load('pancancer_cnv.mat');
CNV_raw_data = CNV_raw.allthresholded;
CNV_raw_gene = cellstr(CNV_raw.gene_information.VarName1(2:end));
CNV_raw_sample = CNV_raw.sample_name;
% load expr
Expr_raw = load('pancancer_expr.mat');
Expr_raw_data = Expr_raw.pancancer_matrix;
Expr_raw_gene = Expr_raw.pancancer_gene;
Expr_raw_sample = Expr_raw.sample_name;
clear CNV_raw Expr_raw
% check NaNs in Expr
num_NaNs_col = sum(isnan(Expr_raw_data),2);
nan_index = find(num_NaNs_col>nan_threshold*length(Expr_raw_sample)); % check if the proportion of NaNs are larger than the threshold
Expr_raw_data(nan_index,:) = []; % remove genes that have too many NaNs
Expr_raw_gene(nan_index) = [];
% fill the NaNs with mean value
for i = 1:size(Expr_raw_data)
    non_nan_tmp = Expr_raw_data(i,find(~isnan(Expr_raw_data(i,:))));
    Expr_raw_data(i,find(isnan(Expr_raw_data(i,:)))) = mean(non_nan_tmp);
end
% keep the sample name the same across CNV and Expr
CNV_raw_sample = cellfun(@(x) x(1:15),CNV_raw_sample,'un',0);
Expr_raw_sample = cellfun(@(x) x{1},Expr_raw_sample,'un',0);
[common_sample,IA,IB] = intersect(CNV_raw_sample,Expr_raw_sample);
CNV_raw_data = CNV_raw_data(:,IA);
Expr_raw_data = Expr_raw_data(:,IB);
% % handle the gene symbols
% CNV_raw_gene = cellfun(@(x) strsplit(x,'|'),CNV_raw_gene,'un',0);
% CNV_gene = cell(length(CNV_raw_gene),2);
% for i = 1:length(CNV_raw_gene)
%     tmp = CNV_raw_gene{i,:};
%     CNV_gene(i,1) = tmp(1);
%     if length(tmp)==2
%         CNV_gene(i,2) = tmp(2);
%     end
% end
% 
% Expr_raw_gene = cellfun(@(x) strsplit(x{1},'|'),Expr_raw_gene,'un',0);
% Expr_gene(:,1) = cellfun(@(x) x{1},Expr_raw_gene,'un',0);
% Expr_gene(:,2) = cellfun(@(x) x{2},Expr_raw_gene,'un',0);
% % align CNV and Expr genes
% CNV_index = [];
% Expr_index = [];
% for i = 1:size(CNV_gene,1)
%     gene = CNV_gene(i,1);
%     ind = find(strcmp(Expr_gene(:,1),gene));
%     if ~isempty(ind)
%         CNV_index = [CNV_index;i];
%         Expr_index = [Expr_index;ind(1)];
%     end
% end
% 
% CNVmatrix = CNV_raw_data(CNV_index,:);
% CNVgene = CNV_gene(CNV_index,:);
% Exprmatrix = Expr_raw_data(Expr_index,:);
% Exprgene = Expr_gene(Expr_index,:);
% save('CNV_Expr_data.mat','CNVmatrix','Exprmatrix','CNVgene','Exprgene','common_sample');




%% parameters for data processing
%        Exp filter     null thres     poisson degree
% v1:   5% (1.64)     max/10               5
% v2:   5% (1.64)       max/5               5
% v3:   2.5% (2)         max/10             5
% v4:   2.5% (2)         max/5               5

paraset = {2,@(x) 5*max(x)/max(x),'poly7'};

filter_multiplier = paraset{1};
null_thres_func = paraset{2};
poisson_deg = paraset{3};


%%%%% Filter the CNV data
load('CNV_Expr_data.mat')


%%% PLOT START %%%%
amp_gene = 'ERBB2';
del_gene = 'PTEN';
index_amp = find(strcmp(CNVgene,amp_gene));
index_del = find(strcmp(CNVgene,del_gene));

% amp
expr_amp = Exprmatrix(index_amp,:);
cnv_amp = CNVmatrix(index_amp,:);
h1=plot(expr_amp,cnv_amp,'x','MarkerSize',20);hold on
u_amp = mean(expr_amp);
std_amp = std(cnv_amp);
xx = linspace(min(expr_amp),max(expr_amp),100);
h4=plot(xx,normpdf(xx,u_amp,std_amp),'b-')
a = ylim;
h2=plot([u_amp+filter_multiplier*std_amp,u_amp+filter_multiplier*std_amp],[a(1),a(2)],'r--');
% plot([u_amp+2*std_amp,u_amp+2*std_amp],[a(1),a(2)],'r-')
index_amp_reliable = intersect(find(expr_amp>u_amp+filter_multiplier*std_amp),find(cnv_amp==2));
h3=plot(expr_amp(index_amp_reliable),cnv_amp(index_amp_reliable),'rx','MarkerSize',20);
hold off
legend({'All CNVs','Expression Distribution','2\sigma Threshold','Significant Events'},'location','northwest');
xlabel('Expression Level');ylabel('Copy Number State')
boldify1
set(gca,'fontsize',18);axis square;set(gca,'ytick',[-2,-1,0,1,2]);h1.LineWidth = 1;
h2.LineWidth = 1;h3.LineWidth = 1;h4.LineWidth = 3;
keyboard
savefig(['./amp_event.fig'])
print(['./amp_event.eps'],'-depsc')
close all

% del
expr_del = Exprmatrix(index_del,:);
cnv_del = CNVmatrix(index_del,:);
h1=plot(expr_del,cnv_del,'x','MarkerSize',20);hold on
u_del = mean(expr_del);
std_del = std(cnv_del);
xx = linspace(min(expr_del),max(expr_del),100);
h4=plot(xx,normpdf(xx,u_del,std_del),'b-');
a = ylim;
h2=plot([u_del-filter_multiplier*std_del,u_del-filter_multiplier*std_del],[a(1),a(2)],'r--')
index_del_reliable = intersect(find(expr_del<u_del-filter_multiplier*std_del),find(cnv_del==-2));
h3=plot(expr_del(index_del_reliable),cnv_del(index_del_reliable),'rx','MarkerSize',20)
hold off
legend({'All CNVs','Expression Distribution','2\sigma Threshold','Significant Events'},'location','best');
xlabel('Expression Level');ylabel('Copy Number State')
boldify1
set(gca,'fontsize',18);axis square;set(gca,'ytick',[-2,-1,0,1,2]);h1.LineWidth = 1;
h2.LineWidth = 1;h3.LineWidth = 1;h4.LineWidth = 3;
keyboard
savefig(['./del_event.fig'])
print(['./del_event.eps'],'-depsc')
close all
keyboard
%%%% PLOT END %%%%
CNVamp = zeros(size(CNVmatrix));
CNVdel = zeros(size(CNVmatrix));
for i = 1:size(CNVmatrix,1)
    med = mean(Exprmatrix(i,:));
    offset = std(Exprmatrix(i,:));
    thres_amp = med + filter_multiplier*offset;
    thres_del = med - filter_multiplier*offset;
    ind_amp = intersect(find(Exprmatrix(i,:)>thres_amp),find(CNVmatrix(i,:)==2));
    ind_del = intersect(find(Exprmatrix(i,:)<thres_del),find(CNVmatrix(i,:)==-2));
    CNVamp(i,ind_amp) = 1;
    CNVdel(i,ind_del) = 1;
end

%% Building Null model
% negative binormal
% Amplification
ob_stat_amp = sum(CNVamp,2);
t_amp = floor(null_thres_func(ob_stat_amp));%quantile(ob_stat_amp,0.9);

pf_truncnbin_amp = @(x,r,p) nbinpdf(x,r,p) ./ nbincdf(t_amp,r,p);
start = [0.5,0.2];
parmhat_amp = mle(ob_stat_amp(find(ob_stat_amp<=t_amp)),'pdf',...
    pf_truncnbin_amp,'start',start,'lower',[0,0]);
theta_amp = length(ob_stat_amp(find(ob_stat_amp<=t_amp)))/length(ob_stat_amp);
p0_amp = theta_amp/nbincdf(t_amp,parmhat_amp(1),parmhat_amp(2));
%parmhat_amp = nbinfit(ob_stat_amp(find(ob_stat_amp<=10)));
pvalue_amp = nbincdf(ob_stat_amp-1e-6,parmhat_amp(1),parmhat_amp(2),'upper');
qvalue_amp = mafdr(pvalue_amp,'BH',true);
x = 0:max(ob_stat_amp);s=p0_amp*nbinpdf(x,parmhat_amp(1),parmhat_amp(2));
% figure;plot(x,s,'-','linewidth',5);hold on
% histogram(ob_stat_amp,'BinMethod','integers','Normalization','probability');
% hold off
% figure;hist(pvalue_amp,1000);
[Y,~] = histcounts(ob_stat_amp,'BinMethod','integers','Normalization','probability');
X = (0:max(ob_stat_amp))';
mdl = fitglm(X,Y,poisson_deg,'Distribution','Poisson');   
yp = predict(mdl,X);
% figure;plot(X,Y,'*',X,yp,'s-',X,s,'r-')
lfdr_for_x_amp = exp(log(s)-log(yp'));
lfdr_amp = zeros(size(ob_stat_amp));
for i = 1:length(ob_stat_amp)
    lfdr_amp(i) = lfdr_for_x_amp(ob_stat_amp(i)+1);
    if lfdr_amp(i)>1
        lfdr_amp(i) = 1;
    end
end
%figure;plot(sort(lfdr_amp))
amp_number = length(find(lfdr_amp<=0.1))



%Deletion
ob_stat_del = sum(CNVdel,2);
t_del = floor(null_thres_func(ob_stat_del));%quantile(ob_stat_del,0.975);
pf_truncnbin_del = @(x,r,p) nbinpdf(x,r,p) ./ nbincdf(t_del,r,p);
start = [1,0.1];
parmhat_del = mle(ob_stat_del(find(ob_stat_del<=t_del)),'pdf',...
    pf_truncnbin_del,'start',start,'lower',[0,0]);
theta_del = length(ob_stat_del(find(ob_stat_del<=t_del)))/length(ob_stat_del);
p0_del = theta_del/nbincdf(t_del,parmhat_del(1),parmhat_del(2));
%parmhat_del = nbinfit(ob_stat_del(find(ob_stat_del<=5)));
pvalue_del = nbincdf(ob_stat_del-1e-6,parmhat_del(1),parmhat_del(2),'upper');
qvalue_del = mafdr(pvalue_del,'BH',true);
x = 0:max(ob_stat_del);s=p0_del*nbinpdf(x,parmhat_del(1),parmhat_del(2));
% figure;plot(x,s,'-','linewidth',5);hold on
% h = histogram(ob_stat_del,'BinMethod','integers','Normalization','probability');
% hold off
% figure;hist(pvalue_del,1000);
[Y,~] = histcounts(ob_stat_del,'BinMethod','integers','Normalization','probability');
X = (0:max(ob_stat_del))';
mdl = fitglm(X,Y,poisson_deg,'Distribution','Poisson');   
yp = predict(mdl,X);
% figure;plot(X,Y,'*',X,yp,'s-',X,s,'r-')
lfdr_for_x_del = s./yp';
lfdr_del = zeros(size(ob_stat_del));
for i = 1:length(ob_stat_del)
    lfdr_del(i) = lfdr_for_x_del(ob_stat_del(i)+1);
    if lfdr_del(i)>1
        lfdr_del(i) = 1;
    end
end
% figure;plot(sort(lfdr_del))
del_number = length(find(lfdr_del<=0.1))



%% 
%%%%% Merge with Mutation based scores
BRCA_mutation = load('./gdac.broadinstitute.org_BRCA-TP.MutSigNozzleReport2CV.Level_4.2016012800.0.0/BRCA_mutation.mat');
GeneMutation = cellstr(BRCA_mutation.gene);
p_mutation = BRCA_mutation.p;
q_mutation  =BRCA_mutation.q;
BRCA_mutation_fdr = load('BRCA_mutation_locfdr.mat');
Gene_lfdr_mutation = BRCA_mutation_fdr.Gene;
lfdr_mutation = BRCA_mutation_fdr.locfdr;

fullGene = union(CNVgene(:,1),GeneMutation);
[~,full_CNV,ind_CNV] = intersect(fullGene,CNVgene(:,1));
[~,full_mut,ind_mut] = intersect(fullGene,GeneMutation(:,1));

p_full_CNV_amp = ones(length(fullGene),1);
p_full_CNV_del = ones(length(fullGene),1);
p_full_mut = ones(length(fullGene),1);
p_full_CNV_amp(full_CNV) = pvalue_amp(ind_CNV);
p_full_CNV_del(full_CNV) = pvalue_del(ind_CNV);
p_full_mut(full_mut) = p_mutation(ind_mut);

q_full_CNV_amp = ones(length(fullGene),1);
q_full_CNV_del = ones(length(fullGene),1);
q_full_mut = ones(length(fullGene),1);
q_full_CNV_amp(full_CNV) = qvalue_amp(ind_CNV);
q_full_CNV_del(full_CNV) = qvalue_del(ind_CNV);
q_full_mut(full_mut) = q_mutation(ind_mut);

lfdr_full_CNV_amp = ones(length(fullGene),1);
lfdr_full_CNV_del = ones(length(fullGene),1);
lfdr_full_mut = ones(length(fullGene),1);
lfdr_full_CNV_amp(full_CNV) = lfdr_amp(ind_CNV);
lfdr_full_CNV_del(full_CNV) = lfdr_del(ind_CNV);
lfdr_full_mut(full_mut) = lfdr_mutation(ind_mut);

tmp = [lfdr_full_CNV_amp,lfdr_full_CNV_del,lfdr_full_mut];
merge_lfdr = min(tmp,[],2);

full_res = table(fullGene,p_full_CNV_amp,p_full_CNV_del,p_full_mut,...
    q_full_CNV_amp,q_full_CNV_del,q_full_mut,...
    lfdr_full_CNV_amp,lfdr_full_CNV_del,lfdr_full_mut,merge_lfdr,...
    'VariableNames',{'Gene','p_amp','p_del','p_mut','q_amp','q_del','q_mut','lfdr_amp',...
    'lfdr_del','lfdr_mut','merged_lfdr'});
full_res = sortrows(full_res,{'merged_lfdr'},'ascend');

mut_res = table(GeneMutation,p_mutation,q_mutation,lfdr_mutation,'VariableNames',...
    {'Gene','p','q','lfdr'});
mut_res = sortrows(mut_res,{'p'},'ascend');
%% save results
keyboard
writetable(full_res,['BRCA__result_merged.xlsx']);
amp_res = table(CNVgene(:,1),ob_stat_amp,pvalue_amp,qvalue_amp, lfdr_amp,'VariableNames',...
    {'gene','count','pvalue','qvalue','locfdr'});
amp_res = sortrows(amp_res,'count','descend');

del_res = table(CNVgene(:,1),ob_stat_del,pvalue_del,qvalue_del,lfdr_del,'VariableNames',...
    {'gene','count','pvalue','qvalue','locfdr'});
del_res = sortrows(del_res,'count','descend');

save(['BRCA_CNV_amp.mat'],'amp_res')
save(['BRCA_CNV_del.mat'],'del_res')
save(['BRCA_mutation.mat'],'mut_res')


%%
