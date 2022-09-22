%% F_confidence_correlations

%close all
clear
clc

load('alldata.mat')

%% set colours
bb_inv = [8,81,156
        49,130,189 
        107,174,214
        189,215,231]./255;
    
rr = [252,174,145
    251,106,74
    222,45,38
    165,15,21]./255;

br = [bb_inv;rr];

%%

for i_sub = 1:length(DATA)
    resp   = DATA{i_sub}.Resp;
    bayesd = STIM{i_sub}.bayesd';
    pme    = STIM{i_sub}.PME_mu.LLR_mu;
    agglom = STIM{i_sub}.centroidadd.LLR_z;
    
    d_types = {'Bayesd', 'mu_{PME}', 'Agglom'};
    
    resps = unique(resp);
    for i_resp = 1:length(resps)
        f_resp = resp==resps(i_resp);
        Bayesd.mean(i_sub,i_resp) = mean(bayesd(f_resp));
        Bayesd.sem(i_sub,i_resp) = std(bayesd(f_resp))./sqrt(sum(f_resp));
        
        PME.mean(i_sub,i_resp) = mean(pme(f_resp));
        PME.sem(i_sub,i_resp) = std(pme(f_resp))./sqrt(sum(f_resp));
        
        Agglom.mean(i_sub,i_resp) = mean(agglom(f_resp));
        Agglom.sem(i_sub,i_resp) = std(agglom(f_resp))./sqrt(sum(f_resp));
    end
end

%%

figd
nsubj = length(STIM);

subplot(1,3,1)
for i_resp = 1:length(resps)
    %e = errorbar(repmat(resps(i_resp),[length(STIM),1]), Bayesd.mean(:,i_resp), Bayesd.sem(:,i_resp)); hold on% [], br(i_resp,:)); hold on
    e = errorbar(resps(i_resp), mean(Bayesd.mean(:,i_resp)), std(Bayesd.mean(:,i_resp))./sqrt(nsubj)); hold on
        e.Color = br(i_resp,:); e.CapSize = 20; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
        xticks([-4:2:4])
        title('Bayes'); ylabel('LLR')
        box off
end

subplot(1,3,2)
for i_resp = 1:length(resps)
    %e = errorbar(repmat(resps(i_resp),[length(STIM),1]), Bayesd.mean(:,i_resp), Bayesd.sem(:,i_resp)); hold on% [], br(i_resp,:)); hold on
    e = errorbar(resps(i_resp), mean(PME.mean(:,i_resp)), std(PME.mean(:,i_resp))./sqrt(nsubj)); hold on
        e.Color = br(i_resp,:); e.CapSize = 20; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
        xticks([-4:2:4])
        title('mu_{PME}'); xlabel('Confidence Rating'); 
        box off
end

subplot(1,3,3)
for i_resp = 1:length(resps)
    %e = errorbar(repmat(resps(i_resp),[length(STIM),1]), Bayesd.mean(:,i_resp), Bayesd.sem(:,i_resp)); hold on% [], br(i_resp,:)); hold on
    e = errorbar(resps(i_resp), mean(Agglom.mean(:,i_resp)), std(Agglom.mean(:,i_resp))./sqrt(nsubj)); hold on
        e.Color = br(i_resp,:); e.CapSize = 20; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
        xticks([-4:2:4])
        title({'Agglomerative', 'Clustering'})
        box off
end

set(gcf, 'position', [100 100 1200 400])
path_figures = '/Users/jennlauralee/Google Drive/WedMock/Causal Inference/Analysis/results/';
saveas(gcf, [path_figures 'conf_corr.png']);