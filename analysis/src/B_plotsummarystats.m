%% calculate and save 'allsummstats'
function B_plotsummarystats(STIM, DATA)

nsubj = length(DATA);

for subjidx = 1:nsubj
    stim = STIM{subjidx};
    data = DATA{subjidx};
    
    resp_feeder = data.Resp_feeder;
    resp_conf   = data.Resp;
    N           = stim.N;
    n1          = stim.N1;
    n0          = stim.N0;
    H           = data.Resp_feeder & stim.Feeder;
    F           = data.Resp_feeder & ~stim.Feeder;
    feeder      = stim.Feeder;

    n1s = unique(n1);
    n0s = unique(n0);
    Ns  = unique(N);
    
    %% respond feeder ns
    % respond feeder, as function of n1
    for i_n1 = 1:length(n1s)
        filter_n1 = n1 == n1s(i_n1);
        npigs{subjidx}.n1(i_n1) = sum(resp_feeder(filter_n1))./length(resp_feeder(filter_n1));
    end

    % respond feeder, as function of n0
    for i_n0 = 1:length(n0s)
        filter_n0 = n0 == n0s(i_n0);
        npigs{subjidx}.n0_C1(i_n0) = sum(resp_feeder(filter_n0 & feeder))./length(resp_feeder(filter_n0 & feeder));
    end

    n0C0s = unique(n0(~feeder));
    for i_n0 = 1:length(n0C0s)
        filter_n0 = n0 == n0C0s(i_n0);
        npigs{subjidx}.n0_C0(i_n0) = sum(resp_feeder(filter_n0 & ~feeder))./length(resp_feeder(filter_n0 & ~feeder));
    end

    % respond feeder, as function of both n1 and n0
    for i_n1 = 1:length(n1s)
        for i_n0 = 1:length(n0s)
            filter_n0n1 = n1 == n1s(i_n1) & n0 == n0s(i_n0);
            npigs{subjidx}.n0n1(i_n0, i_n1) = sum(resp_feeder(filter_n0n1))./length(resp_feeder(filter_n0n1));
        end
    end
    
    % respond feeder, as function of both N and n0, conditioned on C1
    for i_N = 1:length(Ns)
        for i_n0 = 1:length(n0s)
            filter_Nn0 = N == Ns(i_N) & n0 == n0s(i_n0);
            npigs{subjidx}.Nn0_C1(i_N, i_n0) = sum(resp_feeder(filter_Nn0 & feeder))./length(resp_feeder(filter_Nn0 & feeder));
        end
    end
    
    % respond feeder, as function of both N and n1, conditioned on C1
    for i_N = 1:length(Ns)
        for i_n1 = 1:length(n1s)
            filter_Nn1 = N == Ns(i_N) & n1 == n1s(i_n1);
            npigs{subjidx}.Nn1_C1(i_N, i_n1) = sum(resp_feeder(filter_Nn1 & feeder))./length(resp_feeder(filter_Nn1 & feeder));
        end
    end
    
    % respond feeder, as function of N, conditioned on C=0 or C=1
    for i_N = 1:length(Ns)
        filter_N = N == Ns(i_N);
        npigs{subjidx}.N_C1(i_N) = sum(resp_feeder(filter_N & feeder))./length(resp_feeder(filter_N & feeder));
        npigs{subjidx}.N_C0(i_N) = sum(resp_feeder(filter_N & ~feeder))./length(resp_feeder(filter_N & ~feeder));
    end
    
    npigs{subjidx}.Ns  = Ns;
    npigs{subjidx}.n0s = n0s;
    npigs{subjidx}.n1s = n1s;

    %% accuracy

    for i_n = 1:length(Ns)
        hit(i_n) = sum(H(N==Ns(i_n)))./sum(feeder(N==Ns(i_n)));
        false_alarm(i_n) = sum(F(N==Ns(i_n)))./sum(feeder(N==Ns(i_n)));
    end

    accuracy{subjidx}.ns = Ns;
    accuracy{subjidx}.hit = hit;
    accuracy{subjidx}.false_alarm = false_alarm;
    
    accuracy{subjidx}.pcnt  = (sum(resp_feeder(find(feeder))) + sum(~resp_feeder(find(~feeder))))./length(feeder);

    %% density (mean distance)
    
    histedges = [2:1:13];
    quantiles = [0:10:100];
    
    dens{subjidx}        = quantile_presp_conf(stim.MeanDist, quantiles, histedges, STIM, DATA, subjidx);
    %% RF min distance
        
    histedges = [0:0.5:5];
    quantiles = [0:10:100];
    
    mindist{subjidx}     = quantile_presp_conf(stim.MinDist, quantiles, histedges, STIM, DATA, subjidx);
            
    %% RF max local density (gaussian convolve)
   
    histedges = [0:0.5:6];
    quantiles = [0:10:100];
    
    localdens{subjidx}    = quantile_presp_conf(stim.localgaussgrid.max, quantiles, histedges, STIM, DATA, subjidx);

    %% RF mean NN distance
        
    histedges = [0:0.5:10];
    quantiles = [0:10:100];
    
    meanNNdist{subjidx}   = quantile_presp_conf(stim.meanNNdist, quantiles, histedges, STIM, DATA, subjidx);
end

%% Calculate group statistics 

%% RF ns
Npigs = align_ns(npigs);

%% Accuracy
for subjidx = 1:nsubj
    Accuracy.ns(subjidx,:)            = accuracy{subjidx}.ns;
    Accuracy.hit(subjidx,:)           = accuracy{subjidx}.hit;
    Accuracy.false_alarm(subjidx,:)   = accuracy{subjidx}.false_alarm;
    
    
    Accuracy.pcnt(subjidx)            = accuracy{subjidx}.pcnt;
end

Accuracy.hit_mean = mean(Accuracy.hit);
Accuracy.hit_SEM  = std(Accuracy.hit)./sqrt(nsubj);
Accuracy.false_alarm_mean = mean(Accuracy.false_alarm);
Accuracy.false_alarm_SEM  = std(Accuracy.false_alarm)./sqrt(nsubj);

%% RF Density
Dens = calc_quantile_meanSEM(dens);

%% RF Min dist
Mindist = calc_quantile_meanSEM(mindist);

%% RF Local density gaussian
Localdens = calc_quantile_meanSEM(localdens);

%% RF Local density gaussian
MeanNNdist = calc_quantile_meanSEM(meanNNdist);

%% Plots (group)
save(['allsummstats.mat'], 'Npigs', 'Accuracy', 'Dens', 'Mindist', 'Localdens', 'MeanNNdist')
%% RF npigs
Ns = [0:15];

bb15 = [189 215 231;
        176 206 226;
        164 197 221;
        152 188 216;
        140 179 211;
        128 170 206;
        116, 161, 201;
        104, 152, 196;
        92, 143, 191;
        80, 134, 186;
        68, 125, 181;
        56, 116, 176;
        44, 107, 171;
        32, 98, 166;
        20, 89, 161;
        8, 81, 156]./255;

rr15 = [252, 174, 145;
        246, 163, 136;
        240, 152, 128;
        234, 142, 120;
        228, 131, 111;
        223, 121, 103;
        217, 110, 95;
        211, 99, 87;
        205, 89, 78;
        199, 78, 70;
        194, 68, 62;
        188, 57, 54;
        182, 46, 45;
        176, 36, 37;
        170, 25, 29;
        165, 15, 20]./255;

figd
subplot(2,2,1)
plot([0 max(Ns)], [0.5 0.5], 'k--'); hold on
errorbar(Ns, Npigs.n1_mean, Npigs.n1_SEM, 'LineWidth', 2, 'Color', 'k');
    xlim([0 max(Ns)])
    ylim([0 1])
xlabel('n_a_f_f_i_l_i_a_t_e_d')

subplot(2,2,2)
errorbar(n0C0s, Npigs.n0_C0_mean, Npigs.n0_C0_SEM, 'LineWidth', 2, 'Color', [8, 81, 156]./255); hold on
errorbar(Ns, Npigs.n0_C1_mean, Npigs.n0_C1_SEM, 'LineWidth', 2, 'Color', [165, 15, 20]./255); hold on
%plot(n0s, rf_n0,'Color', 'k', 'LineWidth', 4); hold on
plot([0 max(Ns)], [0.5 0.5], 'k--');
legend({'C=0' 'C=1'})
    xlim([0 max(Ns)])
    ylim([0 1])
xlabel('n_u_n_a_f_f_i_l_i_a_t_e_d')

subplot(2,2,3)
for i_n0 = 1:length(Ns)
    idxs = ~isnan(Npigs.n0n1_mean(i_n0,:));
    %errorbar(ns(idxs),RF.n0n1_mean(i_n0,idxs),RF.n0n1_SEM(i_n0,idxs),'LineWidth', 2, 'Color', cm(i_n0,:)); hold on
    plot(Ns(idxs),Npigs.n0n1_mean(i_n0,idxs),'o-','LineWidth', 2, 'Color', bb15(i_n0,:)); hold on
end
xlim([0 max(Ns)])
    ylim([0 1])
    yticks([0 0.5 1])
    colormap(bb15);
cb = colorbar('Ticks', linspace(0,2,length(Ns)),...
         'TickLabels',strsplit(num2str(0:2:14)));
    xlabel('n_a_f_f_i_l_i_a_t_e_d')
    ylabel(cb,'n_u_n_a_f_f_i_l_i_a_t_e_d', 'FontSize', 20)

subplot(2,2,4)
for i_n1 = 1:length(Ns)
    idxs = ~isnan(Npigs.n0n1_mean(:,i_n1));
    %errorbar(ns(idxs),RF.n0n1_mean(idxs,i_n1),RF.n0n1_SEM(idxs,i_n1),'LineWidth', 2, 'Color', cm(i_n1,:)); hold on
    plot(Ns(idxs),Npigs.n0n1_mean(idxs,i_n1),'o-','LineWidth', 2, 'Color', rr15(i_n1,:)); hold on
end
   xlim([0 max(Ns)])
    ylim([0 1])
    yticks([0 0.5 1])
    colormap(bb15);
cb = colorbar('Ticks', linspace(0,2,length(Ns)),...
         'TickLabels',strsplit(num2str(0:2:14)));
    xlabel('n_u_n_a_f_f_i_l_i_a_t_e_d')
    ylabel(cb,'n_a_f_f_i_l_i_a_t_e_d', 'FontSize', 20)
set(gcf, 'Position',  [0, 0, 850, 600])



% Distance plots
Ns = [1:5];
     
bb = [189,215,231
        107,174,214
        49,130,189
        8,81,156]./255;
    
rr = [252,174,145
    251,106,74
    222,45,38
    165,15,21]./255;

fig_dist = figd;
subplot(2,2,1)
plot_datadist('Mean distance', Dens, bb, rr, [2 13]); 

subplot(2,2,2)
plot_datadist('Local density', Localdens, bb, rr, [0 6]);

subplot(2,2,3)
plot_datadist('Mean nearest neighbour distance', MeanNNdist, bb, rr, [0 8]);

subplot(2,2,4)
plot_datadist('Minimum distance', Mindist, bb, rr, [0 6]);


% filename = [path_results 'group/RF.png'];
% saveas(gcf,filename)

Ns = [6 9 12 15]
figd

type_names = {'Mean distance between points', 'Minimum distance between points', 'Local density metric', 'Mean nearest neighbour dist'};
dist_types = {'Dens{i_N}', 'Mindist{i_N}', 'Localdens{i_N}', 'MeanNNdist{i_N}'};
xLims = {[2 13] [0 5] [0 6] [0 10]};
for i_type = 1:4
    subplot(2,2,i_type)
for i_N = 1:length(Ns)
    %% RF Density
    plot_histogram(eval(dist_types{i_type}), type_names{i_type}, xLims{i_type}, rr(i_N,:), bb(i_N,:))
end
end



    %% RF MinDist
    plot_histogram(Mindist{i_N}, 'Minimum distance between points', [0 5])

    %% RF Local density
    plot_histogram(Localdens{i_N}, 'Local density metric', [0 60])

    %% RF mean NN dist
    plot_histogram(MeanNNdist{i_N}, 'Mean NN dist', [0 10])

end

% %% Plots (individual)
%     
% for subjidx = 1:nsubj
%     sub_file_ID = subs{subjidx};
%     %% RF
%     ns = [0:15];
% 
%     N = length(ns);
%     R = linspace(1,0,N);  %// Red from 1 to 0
%     B = linspace(0,1,N);  %// Blue from 0 to 1
%     G = zeros(size(R));   %// Green all zero
% 
%     cm = colormap( [R(:), G(:), B(:)] );  %// create colormap
%     figd
%     subplot(2,2,1)
%     plot([0 max(ns)], [0.5 0.5], 'k:'); hold on
%     plot(rf{subjidx}.n1s, rf{subjidx}.n1,'ko-', 'LineWidth', 4);
%         xlim([0 max(ns)])
%         ylim([0 1])
%     xlabel('n_a_f_f_i_l_i_a_t_e_d')
% 
%     subplot(2,2,2)
%     plot(n0C0s, rf{subjidx}.n0_C0, 'bo-', 'LineWidth', 4); hold on
%     plot(rf{subjidx}.n0s, rf{subjidx}.n0_C1, 'ro-', 'LineWidth', 4); hold on
%     %plot(n0s, rf_n0,'Color', 'k', 'LineWidth', 4); hold on
%     plot([0 max(ns)], [0.5 0.5], 'k:');
%     legend({'C=0' 'C=1'})
%         xlim([0 max(ns)])
%         ylim([0 1])
%     xlabel('n_u_n_a_f_f_i_l_i_a_t_e_d')
% 
%     subplot(2,2,3)
%     cm = colormap( [R(:), G(:), B(:)] );  %// create colormap
%     for i_n0 = 1:length(rf{subjidx}.n0s)
%         idxs = ~isnan(rf{subjidx}.n0n1(i_n0,:));
%         plot(rf{subjidx}.n1s(idxs),rf{subjidx}.n0n1(i_n0,idxs),'-o','Color', cm(i_n0,:)); hold on
%     end
%     xlim([0 max(ns)])
%         ylim([0 1])
%         yticks([0 0.5 1])
%     cb = colorbar('Ticks', linspace(0,2,length(ns)),...
%              'TickLabels',strsplit(num2str(0:2:14)));
%         xlabel('n_a_f_f_i_l_i_a_t_e_d')
%         ylabel(cb,'n_u_n_a_f_f_i_l_i_a_t_e_d', 'FontSize', 20)
% 
%     subplot(2,2,4)
%     for i_n1 = 1:length(rf{subjidx}.n1s)
%         idxs = ~isnan(rf{subjidx}.n0n1(:,i_n1));
%         plot(rf{subjidx}.n0s(idxs),rf{subjidx}.n0n1(idxs,i_n1),'-o','Color', cm(i_n1,:)); hold on
%     end
%        xlim([0 max(ns)])
%         ylim([0 1])
%         yticks([0 0.5 1])
%     cb = colorbar('Ticks', linspace(0,2,length(ns)),...
%              'TickLabels',strsplit(num2str(0:2:14)));
%         xlabel('n_u_n_a_f_f_i_l_i_a_t_e_d')
%         ylabel(cb,'n_a_f_f_i_l_i_a_t_e_d', 'FontSize', 20)
%    
%     filename = [path_results 'individuals/' sub_file_ID '_RF.png'];
%     saveas(gcf,filename)
%     
%     %% Accuracy
%     figd
%     plot(accuracy{subjidx}.ns, accuracy{subjidx}.hit, '-go', 'LineWidth', 4); hold on
%     plot(accuracy{subjidx}.ns, accuracy{subjidx}.false_alarm, '-ro', 'LineWidth', 4);
%     ylim([0 1])
%     xlim([6 15])
%     xlabel('number of pigeons')
%     legend({'hit rate', 'false alarm rate'})
%     
%     filename = [path_results 'individuals/' sub_file_ID '_accuracy{subjidx}.png'];
%     saveas(gcf,filename)
%     
%     %% Density 
%     figd
% 
%     subplot(2,1,1)
%     plot(dens{subjidx}.histcenters, dens{subjidx}.Nhist,  'o-', 'color', 'k'); hold on
%     plot(dens{subjidx}.histcenters, dens{subjidx}.Nhist_C1, 'o-', 'color', 'r'); hold on
%     plot(dens{subjidx}.histcenters, dens{subjidx}.Nhist_C0,  'o-', 'color', 'b');
%     legend({'all','C=1','C=0'})
%     ylabel('N trials')
%     xlim([2 13])
% 
%     subplot(2,1,2)
%     plot(dens{subjidx}.quant_midpts, dens{subjidx}.pRespFeeder, 'o-', 'color', 'k'); hold on
%     plot(dens{subjidx}.quant_midpts_C1, dens{subjidx}.pRespFeeder_C1, 'o-', 'color', 'r'); hold on
%     plot(dens{subjidx}.quant_midpts_C0, dens{subjidx}.pRespFeeder_C0, 'o-', 'color', 'b');
%     legend({'all', 'C=1', 'C=0'})
%     xlabel('Mean distance between points')
%     ylabel('Proportion resp "C=1"')
%     xlim([2 13])
% 
%     title('Proportion responding feeder present')
%     
%     filename = [path_results 'individuals/' sub_file_ID '_dens{subjidx}ity.png'];
%     saveas(gcf,filename)
%     
%     close all
% end

