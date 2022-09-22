%% summary statistics

clear
clc

load('alldata.mat')

% Descriptive stats on data
        % Logistic regression: iv = x-axis, dv = binary
            % p(chat = 1) = logistic(b0 + b1*N + b2*C + b3*X) 
            % where x is the stat of interest

            % per subject meaning you do the logistic regression separately
            % for each subject
            % Then, at cross subject level, do a t-test on each beta 
            
            % Expect significant effect of X on p(C=1) and N for some, etc.
            
            % In short CCN abstract, just report (all p<0.001) etc
            
        % Logistic regression for the simple plot
            % logistic(b0 + b1*N + b2*C)

vars = {'meandist' 'mindist' 'meanNNdist' 'maxlocal'};
for i_sub = 1:length(STIM)
    resp    = categorical(DATA{i_sub}.Resp_feeder);
    N       = STIM{i_sub}.N;
    C       = STIM{i_sub}.Feeder;
    meandist    = STIM{i_sub}.MeanDist;
    mindist     = STIM{i_sub}.MinDist;
    meanNNdist  = STIM{i_sub}.meanNNdist;
    maxlocal    = STIM{i_sub}.localgaussgrid.max';
    % Logistic regression of simple plot
    meas_s = [N C];
    b.simple(i_sub,:) = mnrfit(meas_s,resp)';     % p(chat = 1) = logistic(b0 + b1*N + b2*C)
    
    % Logistic regression of complex plot
    for i_var = 1:length(vars)
        meas = [N C eval(vars{i_var})];             % p(chat = 1) = logistic(b0 + b1*N + b2*C + b3*X) 
        b.complex{i_var}(i_sub,:) = mnrfit(meas,resp)';
    end        
end
% T-test across subs
[tt.simple.t, tt.simple.p] = ttest(b.simple);
for i_var = 1:length(vars)
    [tt.complex.t(i_var,:), tt.complex.p(i_var,:)] = ttest(b.complex{i_var});
end
    