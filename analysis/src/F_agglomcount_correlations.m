%% F_agglomcount_correlations

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

%% kcount by resp

for i_sub = 1:length(DATA)
    count_type = STIM{i_sub}.k_count.unseededcentroidadd;

    
    kcount_    = count_type;
    resp   = DATA{i_sub}.Resp;
    
    resps = unique(resp);
    for i_resp = 1:length(resps)
        f_resp = resp==resps(i_resp);
        
        kcount.mean(i_sub, i_resp) = mean(kcount_(f_resp));
        kcount.sem(i_sub, i_resp) = std(kcount_(f_resp))./sqrt(sum(f_resp));
    end
end

figd
nsubj = length(STIM);

for i_resp = 1:length(resps)
    %e = errorbar(repmat(resps(i_resp),[length(STIM),1]), Bayesd.mean(:,i_resp), Bayesd.sem(:,i_resp)); hold on% [], br(i_resp,:)); hold on
    e = errorbar(resps(i_resp), mean(kcount.mean(:,i_resp)), std(kcount.mean(:,i_resp))./sqrt(nsubj)); hold on
        e.Color = br(i_resp,:); e.CapSize = 20; e.LineWidth = 2; e.MarkerSize = 15;
    lp = plot(resps(i_resp), mean(kcount.mean(:,i_resp)));
        lp.Color = 'k';
        xticks([-4:2:4])
        title('Response'); ylabel('Predicted RT (# of iterations)')
end

%% kcount x RT scatterplot
figd

for i_sub = 1:length(DATA)
    subplot(2,5,i_sub)
    count_type = STIM{i_sub}.k_count.unseededcentroidadd;
    kcount = count_type;
    
    RT = DATA{i_sub}.RT;
    
    [rho(i_sub),pval] = corr(kcount, RT, 'Type', 'Spearman');  
    disp(rho(i_sub));
    disp(pval);
    
    sp = scatter(kcount, RT, 'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0.3);
    text(60,8,['\rho = ' num2str(round(rho(i_sub),2))], 'Color', 'red', 'FontSize', 16);
    text(60,5, ['p = ' num2str(pval)], 'Color', 'red', 'FontSize', 16);

    ax = gca;
    offsetAxes(ax);
    
    ylim([0 10])
    xlim([0 120])
    set(gcf, 'Position', [100 100 1100 400])
end




%% RT correlated with binned LLR(s)

nbins = 10;
quantiles = [0:nbins:100];

for i_sub = 1:length(DATA)
    kcount_    = count_type;
    Bayesd = STIM{i_sub}.bayesd';
    MJP    = STIM{i_sub}.brutemax.posterior.LLR_z_mu;
    Agglom = STIM{i_sub}.centroidadd.LLR_z;
    
    d_types = {'Bayesd'};%, 'MJP', 'Agglom'};
    
    for i_d = 1:length(d_types)
        %subplot(length(d_types),1,i_d)
        llr = eval(d_types{i_d});

        quant = prctile(llr, quantiles);
        quant_midpts{i_d}(i_sub,:) = mean([quant(2:end); quant(1:end-1)]);
                
       for i_quantile = 1:length(quant)-1
        if i_quantile == length(quant)-1 % last quantile includes final entry
            bin(:,i_quantile) = llr>=quant(i_quantile) & llr <=quant(i_quantile+1);
        else 
            bin(:,i_quantile) = llr>=quant(i_quantile) & llr <quant(i_quantile+1);
        end
        
        kcount_binned{i_d}(i_sub,i_quantile) = mean(kcount_(bin(:,i_quantile)));
       end
    end
end

titles = {'Bayes', 'Maximum Joint Posterior', 'Agglomerative Clustering'};

figd
for i_d = 1:length(d_types)
    Kcount_Bin.mean{i_d} = mean(kcount_binned{i_d});
    Kcount_Bin.std{i_d} = std(kcount_binned{i_d})./sqrt(nsubj);
    
    subplot(1,length(d_types),i_d)
    e = errorbar(mean(quant_midpts{i_d}), Kcount_Bin.mean{i_d}, Kcount_Bin.std{i_d}); %hold on
    e.Color = 'k'; e.CapSize = 20; e.LineStyle = '-'; e.LineWidth = 2; e.Marker = '.'; e.MarkerSize = 15;
    title(titles{i_d});
    if i_d == 1
        ylabel('Predicted RT (# of iterations)')
    end
    xlabel('LLR')
    xlim([min(mean(quant_midpts{i_d})) max(mean(quant_midpts{i_d}))])
    
    ax = gca;
    offsetAxes(ax);
    set(gcf, 'Position', [100 100 400 400])
    %ylim([0.2 1.2])
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