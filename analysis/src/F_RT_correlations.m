%% F_RT_correlations

close all
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

%% RT by resp

for i_sub = 1:length(DATA)
    RT_    = DATA{i_sub}.RT;
    resp   = DATA{i_sub}.Resp;
    
    resps = unique(resp);
    for i_resp = 1:length(resps)
        f_resp = resp==resps(i_resp);
        
        RT.mean(i_sub, i_resp) = mean(RT_(f_resp));
        RT.sem(i_sub, i_resp) = std(RT_(f_resp))./sqrt(sum(f_resp));
    end
end

figd
nsubj = length(STIM);

for i_resp = 1:length(resps)
    %e = errorbar(repmat(resps(i_resp),[length(STIM),1]), Bayesd.mean(:,i_resp), Bayesd.sem(:,i_resp)); hold on% [], br(i_resp,:)); hold on
    e = errorbar(resps(i_resp), mean(RT.mean(:,i_resp)), std(RT.mean(:,i_resp))./sqrt(nsubj)); hold on
        e.Color = br(i_resp,:); e.CapSize = 20; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
        xticks([-4:2:4])
        title('Response'); ylabel('RT')
end


%% RT correlated with binned LLR(s)

nbins = 10;
quantiles = [0:nbins:100];

for i_sub = 1:length(DATA)
    RT_    = DATA{i_sub}.RT;
    Bayesd = STIM{i_sub}.bayesd';
    PME    = STIM{i_sub}.PME_mu.LLR_mu;
    Agglom = STIM{i_sub}.centroidadd.LLR_z;
    
    d_types = {'Bayesd'};%, 'PME', 'Agglom'};
    
    for i_d = 1:length(d_types)
        %subplot(3,1,i_d)
        llr = eval(d_types{i_d});

        quant = prctile(llr, quantiles);
        quant_midpts{i_d}(i_sub,:) = mean([quant(2:end); quant(1:end-1)]);
                
       for i_quantile = 1:length(quant)-1
        if i_quantile == length(quant)-1 % last quantile includes final entry
            bin(:,i_quantile) = llr>=quant(i_quantile) & llr <=quant(i_quantile+1);
        else 
            bin(:,i_quantile) = llr>=quant(i_quantile) & llr <quant(i_quantile+1);
        end
        
        RT_binned{i_d}(i_sub,i_quantile) = mean(RT_(bin(:,i_quantile)));
       end
    end
end

titles = {'Bayes', 'PME', 'Agglomerative Clustering'};

figd
for i_d = 1:length(d_types)
    RT_Bin.mean{i_d} = mean(RT_binned{i_d});
    RT_Bin.std{i_d} = std(RT_binned{i_d})./sqrt(nsubj);
    
    subplot(1,length(d_types),i_d)
    e = errorbar(mean(quant_midpts{i_d}), RT_Bin.mean{i_d}, RT_Bin.std{i_d}); %hold on
    e.Color = 'k'; e.CapSize = 20; e.LineStyle = '-'; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
    title(titles{i_d});
    if i_d == 1
        ylabel('Reaction Time (s)')
    end
    xlabel('LLR')
    xlim([min(mean(quant_midpts{i_d})) max(mean(quant_midpts{i_d}))])
    ylim([0.2 1.2])
    
   ax = gca;
   offsetAxes(ax);
   set(gcf, 'Position', [100 100 400 400])
end




% subplot(1,3,2)
% for i_resp = 1:length(resps)
%     %e = errorbar(repmat(resps(i_resp),[length(STIM),1]), Bayesd.mean(:,i_resp), Bayesd.sem(:,i_resp)); hold on% [], br(i_resp,:)); hold on
%     e = errorbar(resps(i_resp), mean(MJP.mean(:,i_resp)), std(MJP.mean(:,i_resp))./sqrt(nsubj)); hold on
%         e.Color = br(i_resp,:); e.CapSize = 20; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
%         xticks([-4:2:4])
%         title('MJP'); xlabel('Confidence Rating'); 
% end
% 
% subplot(1,3,3)
% for i_resp = 1:length(resps)
%     %e = errorbar(repmat(resps(i_resp),[length(STIM),1]), Bayesd.mean(:,i_resp), Bayesd.sem(:,i_resp)); hold on% [], br(i_resp,:)); hold on
%     e = errorbar(resps(i_resp), mean(Agglom.mean(:,i_resp)), std(Agglom.mean(:,i_resp))./sqrt(nsubj)); hold on
%         e.Color = br(i_resp,:); e.CapSize = 20; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
%         xticks([-4:2:4])
%         title({'Agglomerative', 'Clustering'})
% end