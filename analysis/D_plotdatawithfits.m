%% D_plotdatawithfits

function [fig_N, fig_Nn_C1, fig_summstats] = D_plotdatawithfits(STIM, DATA, Modelresults, modelidx, path_summstats, path_results, paffmodel)
    load(path_summstats);
    nsubj = length(DATA);

    ns = [0:15];
    Ns = [6 9 12 15];
    
    bb = [189,215,231
            107,174,214
            49,130,189
            8,81,156]./255;

    rr = [252,174,145
        251,106,74
        222,45,38
        165,15,21]./255;

    % get summ stats from model predictions for plotting
    for subjidx = 1:nsubj
        Stim = STIM{subjidx};
        Data = DATA{subjidx};
        modelpred = Modelresults{modelidx}.modelpred(subjidx,:);
        
        resp_feeder = Data.Resp_feeder;
        N           = Stim.N;
        n1          = Stim.N1;
        n0          = Stim.N0;
        H           = Data.Resp_feeder & Stim.Feeder;
        F           = Data.Resp_feeder & ~Stim.Feeder;
        feeder      = Stim.Feeder;

        n1s = unique(n1);
        n0s = unique(n0);

        %% respond feeder
        
        % respond feeder, as function of N
        for i_N = 1:length(Ns)
            filter_N = N == Ns(i_N);
            mp.RF{subjidx}.N_C1(i_N) = mean(modelpred(filter_N & feeder));
            mp.RF{subjidx}.N_C0(i_N) = mean(modelpred(filter_N & ~feeder));
        end
        
        % respond feeder, as function of n1
        for i_n1 = 1:length(n1s)
            filter_n1 = n1 == n1s(i_n1);
            mp.RF{subjidx}.n1(i_n1) = mean(modelpred(filter_n1));
        end

        % respond feeder, as function of n0
        for i_n0 = 1:length(n0s)
            filter_n0 = n0 == n0s(i_n0);
            mp.RF{subjidx}.n0_C1(i_n0) = mean(modelpred(filter_n0 & feeder));
        end

        n0C0s = unique(n0(~feeder));
        for i_n0 = 1:length(n0C0s)
            filter_n0 = n0 == n0C0s(i_n0);
            mp.RF{subjidx}.n0_C0(i_n0) = mean(modelpred(filter_n0 & ~feeder));
        end

        % respond feeder, as function of both n1 and n0
        for i_n1 = 1:length(n1s)
            for i_n0 = 1:length(n0s)
                filter_n0n1 = n1 == n1s(i_n1) & n0 == n0s(i_n0);
                mp.RF{subjidx}.n0n1(i_n0, i_n1) = mean(modelpred(filter_n0n1 & feeder));
            end
        end
        
        % respond feeder, as function of N and n0, conditioned on C=1
        for i_N = 1:length(Ns)
            for i_n0 = 1:length(n0s)
                filter_Nn0 = N == Ns(i_N) & n0 == n0s(i_n0);
                mp.RF{subjidx}.Nn0_C1(i_N,i_n0) = mean(modelpred(filter_Nn0 & feeder));
            end
        end
        
        % respond feeder, as function of N and n1, conditioned on C=1
        for i_N = 1:length(Ns)
            for i_n1 = 1:length(n1s)
                filter_Nn1 = N == Ns(i_N) & n1 == n1s(i_n1);
                mp.RF{subjidx}.Nn1_C1(i_N,i_n1) = mean(modelpred(filter_Nn1 & feeder));
            end
        end

        mp.RF{subjidx}.n0s = n0s;
        mp.RF{subjidx}.n1s = n1s;
        mp.RF{subjidx}.Ns  = Ns;
        
    end
    % RF
    MP.RF = align_ns(mp.RF);
        
    %% N
    fig_N = figd
    
    HC0 = errorbar(Ns, Npigs.N_C0_mean, Npigs.N_C0_SEM, 'LineWidth', 2, 'Color', bb(4,:)); hold on    
    myshadedarea(Ns, MP.RF.N_C0_mean, MP.RF.N_C0_SEM, 'blue'); hold on
    
    HC1 = errorbar(Ns, Npigs.N_C1_mean, Npigs.N_C1_SEM, 'LineWidth', 2, 'Color', rr(4,:)); hold on    
    myshadedarea(Ns, MP.RF.N_C1_mean, MP.RF.N_C1_SEM, 'red'); hold on
    
    plot([0 max(ns)], [0.5 0.5], 'k--'); hold on
        
    xlim([5 16]);
    xticks(Ns);
    
    ylabel('Proportion resp "C=1"')
    xlabel('N pigeons')
    
    %legend([HC0, HC1], 'C=0', 'C=1')

    set(fig_N, 'Position',  [100, 100, 400, 400])
    
    %% n1 and n0 by N _ C1
    fig_Nn_C1 = figd
    subplot(1,2,1)
    plot([0 max(ns)], [0.5 0.5], 'k:'); hold on
    
    for i_N = 1:length(Ns)
        indxs{i_N} = find(~isnan(mean(squeeze(Npigs.Nn1_C1(:,i_N,:)))));
        HC1(i_N) = errorbar(ns(indxs{i_N}), Npigs.Nn1_C1_mean(i_N,indxs{i_N}), Npigs.Nn1_C1_SEM(i_N,indxs{i_N}), 'LineStyle', 'none', 'LineWidth', 2, 'Color', rr(i_N,:)); hold on
        
        myshadedarea(ns(indxs{i_N}), MP.RF.Nn1_C1_mean(i_N,indxs{i_N}), MP.RF.Nn1_C1_SEM(i_N,indxs{i_N}), rr(i_N,:)); hold on
    end
        xlim([0 max(ns)])
        ylim([0 1])
    xlabel('n_a_f_f_i_l_i_a_t_e_d')
    %legend([HC1], num2str(Ns'))
    ylabel('Proportion resp "C=1"')
    title('C=1')
    
    subplot(1,2,2)
    plot([0 max(ns)], [0.5 0.5], 'k:'); hold on
    
    for i_N = 1:length(Ns)
        indxs{i_N} = find(~isnan(mean(squeeze(Npigs.Nn0_C1(:,i_N,:)))));
        HC1(i_N) = errorbar(ns(indxs{i_N}), Npigs.Nn0_C1_mean(i_N,indxs{i_N}), Npigs.Nn0_C1_SEM(i_N,indxs{i_N}), 'LineStyle', 'none', 'LineWidth', 2, 'Color', bb(i_N,:)); hold on
        
        myshadedarea(ns(indxs{i_N}), MP.RF.Nn0_C1_mean(i_N,indxs{i_N}), MP.RF.Nn0_C1_SEM(i_N,indxs{i_N}), bb(i_N,:)); hold on
    end
        xlim([0 max(ns)])
        ylim([0 1])
    xlabel('n_u_n_a_f_f_i_l_i_a_t_e_d')
    title('C=1')
    %legend([HC1], num2str(Ns'))
    set(fig_Nn_C1, 'Position',  [100, 100, 875, 400])
    

    %% distance-related heuristics, MP

    MP.Dens = quantile_MP_presp_conf(Dens, Modelresults, modelidx, STIM);
    MP.Localdens = quantile_MP_presp_conf(Localdens, Modelresults, modelidx, STIM);
    MP.MeanNNdist =  quantile_MP_presp_conf(MeanNNdist, Modelresults, modelidx, STIM);
    MP.Mindist    =  quantile_MP_presp_conf(Mindist, Modelresults, modelidx, STIM);
    
    fig_summstats = figd;
    
    plot_MPdist('Mean distance', Dens, MP.Dens, bb, rr, [3 12.5],[4,1],1);
    plot_MPdist('Local density', Localdens, MP.Localdens, bb, rr, [0 4.5],[4,1],2);
    plot_MPdist('Mean nearest neighbour distance', MeanNNdist, MP.MeanNNdist, bb, rr, [1 7.5],[4,1],3);
    plot_MPdist('Minimum distance', Mindist, MP.Mindist, bb, rr, [0 5.5],[4,1],4);
    
    set(gcf,'position',[100 100 1250 275])
    % Append MP (model predictions) to Modelresults
    Modelresults{modelidx}.MP = MP;
    if paffmodel
        save([path_results 'results_' num2str(modelidx) '_paff.mat'], 'Modelresults', 'DATA', 'STIM');
    else
        save([path_results 'results_' num2str(modelidx) '.mat'], 'Modelresults', 'DATA', 'STIM');
    end
end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     %% n
%     fig_n = figd
%     subplot(1,2,1)
%     plot([0 max(ns)], [0.5 0.5], 'k:'); hold on
%     errorbar(ns, Npigs.n1_mean, Npigs.n1_SEM, 'LineWidth', 2, 'Color', 'r'); hold on
% 
%     indxs =  find(~isnan(MP.RF.n1_mean));
%     myshadedarea(ns(indxs), MP.RF.n1_mean(indxs), MP.RF.n1_SEM(indxs), 'red');
%         xlim([0 max(ns)])
%         ylim([0 1])
%     xlabel('n_a_f_f_i_l_i_a_t_e_d')
% 
%     subplot(1,2,2)
%     HC0 = errorbar(n0C0s, Npigs.n0_C0_mean, Npigs.n0_C0_SEM, 'LineWidth', 2, 'Color', 'b'); hold on
%     myshadedarea(n0C0s', MP.RF.n0_C0_mean, MP.RF.n0_C0_SEM, 'blue'); hold on
% 
%     HC1 = errorbar(ns, Npigs.n0_C1_mean, Npigs.n0_C1_SEM, 'LineWidth', 2, 'Color', 'r'); hold on
%     indxs =  find(~isnan(MP.RF.n0_C1_mean));
%     myshadedarea(ns(indxs), MP.RF.n0_C1_mean(indxs), MP.RF.n0_C1_SEM(indxs), 'red'); hold on
%     plot([0 max(ns)], [0.5 0.5], 'k:');
%     legend([HC0, HC1], 'C=0', 'C=1')
%         xlim([0 max(ns)])
%         ylim([0 1])
%     xlabel('n_u_n_a_f_f_i_l_i_a_t_e_d')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%      
