  
clear
clc
load('results_1_paff.mat')
paffmodel = 1;
modelidx = 1;

for subjidx = 1:10
    par_est = Modelresults{1}.par_est(subjidx, :);
    stim = STIM{subjidx};
    [~, ~, false_d{subjidx}]       = C_modelpredictions_v2(par_est, modelidx, stim, [], paffmodel);
end