%% D_callplotdatawfits

close all
clear
clc

addpath(genpath('/Users/jennlauralee/GitHub Repos/pigeons/Analysis/'));

path_summstats = '/Users/jennlauralee/GitHub Repos/pigeons/Analysis/allsummstats.mat';
path_figures = '/Users/jennlauralee/GitHub Repos/pigeons/Analysis/results/model fits/';
path_results = '/Users/jennlauralee/GitHub Repos/pigeons/Analysis/results/';

load('alldata.mat')

paffmodel = 1;
if paffmodel
    pafftext = '_paff';
    pafftitle = 'P';
else
    pafftext = '';
    pafftitle = '';
end

%bayesAC = [1 8];

% bayesAC = [1 52 39]; 
% bestmodels = [49 40 47];
% % OLDallmodels = fliplr([2,3,4,5,14,16,... %H
% %           41:47,... %D+
% %           39,40,... %C+
% %           48,49,50,51,... %B+
% %           53,52,1]); %A
% 
% allmodels = fliplr([...
%     2,3,4,5,14,16,...           %H
%     43,44,45,46,47,... %D+A   
%     39,40,...%C+A
%     49,50,51,...%B+A
%     53,52,1]); %A

allmodels = [1:19];
allpaffmodels = [1:2,4:13];
bayesagglom = [1,8];
As = [1:3];
Apaff = [1:2];
Bs = [4:6];
Cs = [7,8];
Ds = [9:13];
Hs = [14:19];
    


allmodelnames = {'A1', 'A2', 'A3',...
    'B4', 'B5', 'B6', ...
    'C7', 'C8',...
    'D9', 'D10', 'D11', 'D12', 'D13', ...
    'H14', 'H15', 'H16', 'H17', 'H18', 'H19'};

%% Set which models to plot
%modelgroups = {As};%{As,Bs,Cs,Ds,Hs};
%modelgroups = {Apaff,Bs,Cs,Ds};
modelgroups = {allpaffmodels};%= {As,Bs,Cs,Ds,Hs};
whichfig = 2; % 1 = Ns, 2 = NC1C0, 3 = summstats


for i_modelgroup = 1:length(modelgroups)
    models = modelgroups{i_modelgroup};
    %% fig_Ns
    switch whichfig
    case 1

        figNs = figd;
        subplotdim = [ceil(length(models)/3),6];
        set(gcf,'Position',[100 100 300*subplotdim(2) 200*subplotdim(1)])

        i_plot = 1;
        for i_m = 1:length(models)
            modelidx = models(i_m);

            load(['results_' num2str(modelidx) pafftext '.mat'])

            plot_fig_Ns(subplotdim, i_plot, STIM, DATA, Modelresults, modelidx, allmodelnames, paffmodel, path_summstats, path_results)
            i_plot = i_plot + 2;
            box off
        end

    %     if length(models) == 3
    %         set(gcf, 'position', [100 100 500 600])
    %         saveas(figNs, [path_figures pafftext 'best_fig_N.svg']);
    %     elseif length(models) == 2
    %         set(gcf, 'position', [100 100 500 400])
    %         saveas(figNs, [path_figures pafftext 'main_fig_N.svg']);
    %     end

    %% fig_NC1C0s
    case 2
        
        fig_NC1C0s = figd;
        subplotdim = [ceil(length(models)/5),5];
        set(gcf,'Position',[100 100 600*subplotdim(2) 400*subplotdim(1)])

        for i_m = 1:length(models)
            modelidx = models(i_m);

            load(['results_' num2str(modelidx) pafftext '.mat'])
            plot_fig_N_C1C0(subplotdim, i_m, STIM, DATA, Modelresults, modelidx, path_summstats)
                set(gca,'FontSize', 14)


            title([allmodelnames{modelidx} pafftitle],'FontSize',16);

            box off

        end

    %     if isequal(modelname,'bayesAC')
    %         set(gcf, 'position', [100 100 900 350])
    %         saveas(fig_NC1C0s, [path_figures pafftext 'bayesAC_NC1C0.svg']);
    %     elseif length(models) == 2
    %         set(gcf,'Position', [100 100 600 350])
    %         saveas(fig_NC1C0s, [path_figures pafftext 'main_fig_NC1C0.svg']);
    %     elseif isequal(modelname,'allmodels')
    %         set(gcf,'Position', [100 100 800 1000])
    %         saveas(fig_NC1C0s, [path_figures pafftext 'all_fig_NC1C0.svg']);
    %     end
    %         
    case 3 
        figsummstats = figd;

        if length(models) == 3
            subplotdims = [3,4];
        elseif length(models) == 2
            subplotdims = [2,4];
        elseif length(models) == 1
            subplotdims = [1,4];
        elseif length(models)>3
            error('too many models for summstat plot, max 3');
        end

        i_plot = 1;
        for i_m = 1:length(models)
            modelidx = models(i_m);

            load(['results_' num2str(modelidx) pafftext '.mat'])

            plot_fig_summstats(subplotdims, i_plot, STIM, DATA, Modelresults, modelidx, allmodelnames, paffmodel, path_summstats)        
            i_plot = i_plot + 4;
        end


        set(gcf, 'position', [100 100 1000 300*length(models)])
    end
    
end
    