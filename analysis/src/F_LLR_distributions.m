%% G_LLR_distributions

% For bayes and AC
close all
clear
clc

bb = [189,215,231
        107,174,214
        49,130,189
        8,81,156]./255;

rr = [252,174,145
    251,106,74
    222,45,38
    165,15,21]./255;

% bayesAC = [1 39];
% bestmodels = [49 40 47];

models = [1,100];

load('alldata.mat')

Ns = [6 9 12 15];

for i_m = 1:length(models)
    modelidx = models(i_m);
    for i_N = 1:length(Ns)
        nsubj = length(STIM);
         for subjidx = 1:nsubj
 
            filt_N = STIM{subjidx}.N == Ns(i_N);
            feeder = STIM{subjidx}.Feeder;

            switch modelidx
                case 1
                    llr = STIM{subjidx}.bayesd;
%                 case 39
%                     llr = STIM{subjidx}.centroidadd.LLR_z;
                case 9
                    llr = STIM{subjidx}.brutemax.posterior.LLR_z_mu;
%                 case 40
%                     llr = STIM{subjidx}.mostlikelyz.LLR;
%                 case 49
%                     llr = STIM{subjidx}.maxmargmu.posterior.LLR;
                case 100
                   llr = STIM{subjidx}.bayesd_07paff;
            end


            llr_N = llr(filt_N);
            llr_N_C1 = llr(filt_N & feeder);
            llr_N_C0 = llr(filt_N & ~feeder);

            LLR{modelidx}{i_N}.all(:,subjidx)= llr_N;
            LLR{modelidx}{i_N}.C1(:,subjidx) = llr_N_C1;
            LLR{modelidx}{i_N}.C0(:,subjidx) = llr_N_C0;
         end
         
            LLR{modelidx}{i_N}.all  = LLR{modelidx}{i_N}.all(:);
            LLR{modelidx}{i_N}.C1   = LLR{modelidx}{i_N}.C1(:);
            LLR{modelidx}{i_N}.C0   = LLR{modelidx}{i_N}.C0(:);
    end
end

%% plot
figd

%Nbins = 15;
binedges = [-10:1:25];
for i_m = 1:length(models)
    modelidx = models(i_m);
    LLR_ = LLR{modelidx};
    subplot(length(models),1,i_m)
    for i_N = 1:length(Ns)
        [N_all, edges_all] = histcounts(LLR_{i_N}.all, binedges);
        [N_1, edges_1] = histcounts(LLR_{i_N}.C1, binedges);
        [N_0, edges_0] = histcounts(LLR_{i_N}.C0, binedges);

        centers_all = mean([edges_all(1:end-1); edges_all(2:end)]);
        centers_1 = mean([edges_1(1:end-1); edges_1(2:end)]);
        centers_0 = mean([edges_0(1:end-1); edges_0(2:end)]);

        %plot(centers_all, N_all, 'o-'); hold on
        plot(centers_1, N_1, '-', 'Color', rr(i_N,:)); hold on
        plot(centers_0, N_0, '-', 'Color', bb(i_N,:)); hold on
        plot([0 0], [0 1000], '--', 'Color',' k');
        xlim([-10 25])
    end
end
set(gcf, 'Position', [100 100 700 600])
%set(gcf, 'Position', [100 100 700 800])
