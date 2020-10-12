%%%%% Script for plotting simulation data
% Author: Le Yang
% This script is used to plot generated simulation 
% results for subnetwork detection in the paper:
%   FDRnet: A Novel Method to Identify Signicantly 
%   Perturbed Subnetworks in Cancer

clear all;close all;clc
networks = {'irefindex9','biogrid','reactome'};
beta = [0.11,0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01];


for n = 1:length(networks)
    %% FDRnet
    result_fdrnet = load(['../FDRnet/',networks{n},'fdr_scores.mat']);
    FDRnet_f = result_fdrnet.result_scores_post_neighbor{1}.fscore;
    FDRnet_fsub = result_fdrnet.result_scores_post_neighbor{1}.fsub;
    %% BioNet
    result_bionet = load(['../BioNet/',networks{n},'_bionet_score.mat']);
    BioNet_f = result_bionet.bionet_f;
    BioNet_fsub = result_bionet.bionet_fsub;
    %% hHotNet
    result_hhotnet = load(['../hHotNet/hierarchical-hotnet-master/',networks{n},'_hhotnet_score_v100_full.mat']);
    hHotNet_f = result_hhotnet.hhotnet_f;
    hHotNet_fsub = result_hhotnet.hhotnet_fsub;
    %% HotNet2
    result_hotnet2 = load(['../HotNet2/',networks{n},'_hotnet2_score.mat']);
    hotnet2_f = result_hotnet2.hotnet2_f;
    hotnet2_fsub = result_hotnet2.hotnet2_fsub;
    %% ClustEx
    result_clustex = load(['../ClustEx/',networks{n},'_clustex_score.mat']);
    clustex_f = result_clustex.clustex_f;
    clustex_fsub = result_clustex.clustex_fsub;
    %% RegMOD
    result_regmod = load(['../RegMOD/',networks{n},'_regmod_score.mat']);
    regmod_f = result_regmod.regmod_f;regmod_f(find(isnan(regmod_f))) = 0;
    regmod_fsub = result_regmod.regmod_fsub;
    %% Plot f score: each result is a 11*10 matrix
    marker_s = 10;
    figPos = get(0,'defaultfigureposition');
    width = 660;
    height = 660;
    figure('Position', [figPos(1), figPos(2), width, height]);hold on
    errorbar(beta,mean(FDRnet_f',1),std(FDRnet_f',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
        'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
     errorbar(beta,mean(hotnet2_f',1),std(hotnet2_f',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
         'MarkerEdgeColor',' red','MarkerFaceColor','w', 'linewidth', 2);
     errorbar(beta,mean(clustex_f',1),std(clustex_f',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
         'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
     errorbar(beta,mean(regmod_f',1),std(regmod_f',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
         'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
    errorbar(beta,mean(BioNet_f',1),std(BioNet_f',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
        'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
    errorbar(beta,mean(hHotNet_f',1),std(hHotNet_f',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
        'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
    hold off
    method = {'FDRnet','HotNet2','ClustEx','RegMOD','BioNet','hHotNet'};
    ylim([0,1])
    xlim([0.01,0.11])
    legend(method)
    boldify1;set(gca,'FontSize',18);xticks(sort(beta));box on;axis square
    ylabel('F-Score','FontSize',18)
    xlabel('\alpha','FontSize',18)
    savefig([networks{n},'_f.fig'])
    print(gcf,[networks{n},'_f'],'-depsc')
    close all
    %% Plot fsub score
    marker_s = 10;
    figPos = get(0,'defaultfigureposition');
    width = 660;
    height = 660;
    figure('Position', [figPos(1), figPos(2), width, height]);hold on
    errorbar(beta,mean(FDRnet_fsub',1),std(FDRnet_fsub',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
        'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
     errorbar(beta,mean(hotnet2_fsub',1),std(hotnet2_fsub',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
         'MarkerEdgeColor',' red','MarkerFaceColor','w', 'linewidth', 2);
     errorbar(beta,mean(clustex_fsub',1),std(clustex_fsub',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
         'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
     errorbar(beta,mean(regmod_fsub',1),std(regmod_fsub',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
         'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
    errorbar(beta,mean(BioNet_fsub',1),std(BioNet_fsub',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
        'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
    errorbar(beta,mean(hHotNet_fsub',1),std(hHotNet_fsub',0,1),'-s','CapSize',8, 'MarkerSize',marker_s,...
        'MarkerEdgeColor','red','MarkerFaceColor','w', 'linewidth', 2);
    hold off
    method = {'FDRnet','HotNet2','ClustEx','RegMOD','BioNet','hHotNet'};%{'FDRnet','HotNet2','ClustEx','RegMOD','BioNet','hHotNet'};
    xlim([0.01,0.11])
    legend(method)
    boldify1;set(gca,'FontSize',18);xticks(sort(beta));box on;axis square
    ylabel('Fsub Score','FontSize',18)
    xlabel('\alpha','FontSize',18)
    savefig([networks{n},'_fsub.fig'])
    print(gcf,[networks{n},'_fsub'],'-depsc')
    close all
    %% plot fdr
    %% plot real/estimated FDRs: Figure 3
    for p = 1;%15
        FDRs = [];Group=[];Size = [];%p=12;%a=0.06
        % generate fdr, group and size data
        % group: FDRnet,...
        % size: subnetwork size
        for permute = 1:10
            % FDRnet
            fdrnet_lfdr = result_fdrnet.result_scores_post_neighbor{1}.lfdr{p,permute};
            fdrnet_real_fdr = result_fdrnet.result_scores_post_neighbor{1}.fdr_real{p,permute};
            fdrnet_size = result_fdrnet.result_scores_post_neighbor{1}.size{p,permute};
            ss = length(fdrnet_lfdr);
            FDRs = [FDRs;fdrnet_lfdr];
            tmp_group = -0.2+ones(ss,1);
            Group=[Group;tmp_group];
            Size = [Size;fdrnet_size];
            FDRs = [FDRs;fdrnet_real_fdr];
            tmp_group = 0.2+ones(ss,1);
            Group=[Group;tmp_group];
            Size = [Size;fdrnet_size];
            % BioNet
            bionet_fdr = result_bionet.bionet_fdr;
            bionet_fdr_real = result_bionet.bionet_fdr_real;
            bionet_result = result_bionet.bionet_result.result;
            if ~isempty(bionet_fdr{p,permute})
                FDRs = [FDRs;bionet_fdr{p,permute}];
                tmp_group = -0.2+2*ones(length(bionet_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                tmp_size = cell2mat(cellfun(@(x) length(x),bionet_result{p,permute},'un',0));
                Size = [Size;tmp_size];
                FDRs = [FDRs;bionet_fdr_real{p,permute}];
                tmp_group = 0.2+2*ones(length(bionet_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                Size = [Size;tmp_size];
            end
            % HotNet2
            hotnet_fdr = result_hotnet2.hotnet2_fdr;
            hotnet_fdr_real = result_hotnet2.hotnet2_real_fdr;
            hotnet_result = result_hotnet2.hotnet2_result;
            if ~isempty(hotnet_fdr{p,permute})
                FDRs = [FDRs;hotnet_fdr{p,permute}];
                tmp_group = -0.2+3*ones(length(hotnet_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                tmp_size = cell2mat(cellfun(@(x) length(x),hotnet_result{p,permute},'un',0));
                Size = [Size;tmp_size];
                FDRs = [FDRs;hotnet_fdr_real{p,permute}];
                tmp_group = 0.2+3*ones(length(hotnet_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                Size = [Size;tmp_size];
            end
            % hHotNet
            hhotnet_fdr = result_hhotnet.hhotnet_fdr;
            hhotnet_fdr_real = result_hhotnet.hhotnet_real_fdr;
            hhotnet_result = result_hhotnet.hhotnet_result;
            if ~isempty(hhotnet_fdr{p,permute})
                FDRs = [FDRs;hhotnet_fdr{p,permute}];
                tmp_group = -0.2+4*ones(length(hhotnet_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                tmp_size = cell2mat(cellfun(@(x) length(x),hhotnet_result{p,permute},'un',0));
                Size = [Size;tmp_size];
                FDRs = [FDRs;hhotnet_fdr_real{p,permute}];
                tmp_group = 0.2+4*ones(length(hhotnet_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                Size = [Size;tmp_size];
            end
            % RegMOD
            regmod_fdr = result_regmod.regmod_fdr;
            regmod_fdr_real = result_regmod.regmod_real_fdr;
            regmod_result = result_regmod.regmod_result;
            FDRs = [FDRs;regmod_fdr{p,permute}];
            tmp_group = -0.2+5*ones(length(regmod_fdr{p,permute}),1);
            Group=[Group;tmp_group];
            tmp_size = cell2mat(cellfun(@(x) length(x),result_regmod.regmod_networks{p,permute},'un',0));
            Size = [Size;tmp_size];
            FDRs = [FDRs;regmod_fdr_real{p,permute}];
            tmp_group = 0.2+5*ones(length(regmod_fdr{p,permute}),1);
            Group=[Group;tmp_group];
            Size = [Size;tmp_size];
            % ClustEx
            clustex_fdr = result_clustex.clustex_fdr;
            clustex_fdr_real = result_clustex.clustex_real_fdr;
            clustex_result = result_clustex.clustex_result;
            if ~isempty(clustex_fdr{p,permute})
                FDRs = [FDRs;clustex_fdr{p,permute}];
                tmp_group = -0.2+6*ones(length(clustex_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                tmp_size = cell2mat(cellfun(@(x) length(x),clustex_result{p,permute},'un',0));
                Size = [Size;tmp_size];
                FDRs = [FDRs;clustex_fdr_real{p,permute}];
                tmp_group = 0.2+6*ones(length(clustex_fdr{p,permute}),1);
                Group=[Group;tmp_group];
                Size = [Size;tmp_size];
            end
        end
        
        
        figPos = get(0,'defaultfigureposition');
        width = 920;
        height = 350;
        eFDRs_ind = [find(Group==0.8);find(Group==1.8);find(Group==2.8);find(Group==3.8);find(Group==4.8);find(Group==5.8)];
        rFDRs_ind = [find(Group==1.2);find(Group==2.2);find(Group==3.2);find(Group==4.2);find(Group==5.2);find(Group==6.2)];
        Real = zeros(length(Group),1);Real(eFDRs_ind)=1;Real(rFDRs_ind)=2;
        % size chart: subnetwork size, plot size, plot label
        size_chart = {10,50,'< \fontname{helvetica}  10';...
            50, 100,'< \fontname{helvetica}  50';...
            100,200,'< \fontname{helvetica}  100';...
            100,400,'\geq \fontname{helvetica} 100'};
        Size_plot = zeros(size(Size));
        for i = 1:size(size_chart,1)
            if i == 1
                Size_plot(find(Size<size_chart{i,1})) = size_chart{i,2};
            elseif i == size(size_chart,1)
                Size_plot(find(Size>=size_chart{i,1})) = size_chart{i,2};
            else
                Size_plot(intersect(find(Size<size_chart{i,1}),find(Size>=size_chart{i-1,1}))) = size_chart{i,2};
            end
        end
        
        figure('Position', [figPos(1), figPos(2), width, height]);
        % plot first axes and legend 1
        hAx(1) = axes();hold(hAx(1));hp=[0;0];
        hp(1) = scatter(Group(eFDRs_ind),FDRs(eFDRs_ind),Size_plot(eFDRs_ind),'ro','linewidth',0.001,'jitter','on','jitteramount',0.1,'Parent',hAx(1));
        hp(2) = scatter(Group(rFDRs_ind),FDRs(rFDRs_ind),Size_plot(rFDRs_ind),'bo','linewidth',0.001,'jitter','on','jitteramount',0.1,'Parent',hAx(1));
        line([0.5,6.5],[0.1,0.1],'Linestyle','--','Color','k','Parent',hAx(1))
        %title(['\alpha = ',num2str(beta(p))],'Fontweight','normal')
        legend(hp,{'Estimated FDRs','Exact FDRs'},'location','northwest','FontSize',18,'autoupdate','off')
        ylabel(hAx(1),'False Discovery Rate')
        xlim(hAx(1),[0.5,6.5])
        set(hAx(1),'Xtick',[1,2,3,4,5,6]);
        set(hAx(1),'xticklabel',{'FDRnet','BioNet','HotNet2','hHotNet','RegMOD','ClustEx'})
        box on;boldify1;set(hAx(1),'FontSize',18);hold(hAx(1),'off')
        % plot second empty axes and legend 2
        % a. size chart transformation
        hAx(2) = copyobj(hAx(1),gcf); delete(get(hAx(2),'Children')); hold(hAx(2),'on');
        legend2 = cell(size(size_chart,1),1);
        tmp_hh = zeros(size(size_chart,1),1);
        for i = 1:size(size_chart,1)
            tmp_hh(i) = plot(nan,nan,'ko','markersize',sqrt(size_chart{i,2}),'MarkerFaceColor','w','Parent',hAx(2));
            legend2{i} = size_chart{i,3};
        end
        hold(hAx(2),'off')
        set(hAx(2), 'Color','none', 'XTick',[], ...
            'YAxisLocation','left', 'Box','off','YTick',[])  % copy axes 1 and set it invisible
        hAx(2).YLabel.String='';
        pp=legend(tmp_hh,legend2,'FontSize',18,'location',[0.140 0.501 0.107 0.234],'box','off');
        set(pp,'interpreter','tex')
        % t = table(Group,FDRs,Size,Real);
        % writetable(t,['fdr_plot_data',num2str(beta(p)),'.txt'])
        
        savefig([networks{n},'_fdr_simu',num2str(beta(p)),'notitle_v50.fig'])
        print(gcf,[networks{n},'_fdr_simu',num2str(100*beta(p)),'notitle_v50.eps'],'-depsc','-loose')%,'-opengl')
        close all
end
    %% plot running time
    bionet_time_biogrid = 2097 + 215;
    bionet_time_reactome =  1417+101;
    fdrnet_time_biogrid = 1413.1 + 215;
    fdrnet_time_reactome = 1312 + 101;
    regmod_time_biogrid = 2771 + 84;
    regmod_time_reactome = 1164 + 31;
    clustex_time_biogrid = 36419+3860;
    clustex_time_reactome = 17161 + 2972;
    hotnet2_time_biogrid = (19000-215)*100+63995;
    hotnet2_time_reactome = (10000-101)*100+49089;
    hhotnet_time_biogrid = (213343-5*215)*5 + 20*60;
    hhotnet_time_reactome = (15082-20*101)*20 + 20*60;
    
    time_all = [fdrnet_time_biogrid,fdrnet_time_reactome;...
                bionet_time_biogrid,bionet_time_reactome;...
                hotnet2_time_biogrid,hotnet2_time_reactome;...
                hhotnet_time_biogrid,hhotnet_time_reactome;...
                regmod_time_biogrid,regmod_time_reactome;...
                clustex_time_biogrid,clustex_time_reactome];
    cc= cbrewer('qual','Set1',6);
    b = bar(time_all');
    for i = 1:length(b)
        set(b(i),'FaceColor',cc(i,:));
    end
    method = {'FDRnet','BioNet','HotNet2','hHotNet','RegMOD','ClustEx'};
    legend(method,'location','best','fontsize',16,'position',[0.758,0.633,0.137,0.269]);
    % h=boxplot(time_all');
    set(gca,'Yscale','log')
    % set(h,{'linew'},{1})
    % set(gca,'Xticklabel',{'FDRnet','BioNet','RegMOD','ClustEx','HotNet2','hHotNet'})
    % ylabel('Running time /s')
    % boldify1
    set(gca,'Xtick',[1,2])
    set(gca,'Xticklabel',{"    BioGRID\newline(17501 nodes)","  ReactomeFI\newline(13360 nodes)"})
    ylim([10,5000000])
    ylabel("Running Time (log10(sec))")
    set(gca,'Ytick',[60,600,3600,3600*24,3600*24*2,3600*24*4,3600*24*8,3600*24*16]);grid off;box on
    set(gca,'yticklabel',{'1 min','10 mins','1 hour','1 day', '2 days','4 days','8 days','16 days'});
    xlim([0.5,2.5])
    boldify1
    set(gca,'FontSize',18);
    set(gca,'YGrid','on')
    keyboard
    savefig('speed.fig')
    print(gcf,'speed','-depsc')
            
            
end


