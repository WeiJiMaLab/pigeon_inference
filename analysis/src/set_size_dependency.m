%% set size dependency?

stim = STIM{1};

N = stim.N;
feeder = stim.Feeder

figd
for i_N = 1:4
    N_ = N(i_N);
    subplot(4,1,i_N)
    [count1, edges1] = histcounts(stim.MeanDist_Mu.LLR_z_mu(N == N_ & stim.Feeder));
    [count0, edges0] = histcounts(stim.MeanDist_Mu.LLR_z_mu(N == N_ & ~stim.Feeder));
    
    centers1 = mean([edges1(1:end-1); edges1(2:end)]); 
    centers0 = mean([edges0(1:end-1); edges0(2:end)]);
    
    plot(centers1,count1); hold on
    plot(centers0,count0)
end
    