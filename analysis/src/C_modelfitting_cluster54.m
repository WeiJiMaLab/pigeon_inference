%% C_modelfitting

function C_modelfitting_cluster54(i_job)

[subjidx, run] = find(reshape([1:100],[10,10])== i_job);

[subjidx run]

modelidx = 54;
% Todo: save general settings above subject level in superstructure e.g.
% modelresults{modelidx}. and modelresults{modelidx}.bysubject{subjidx}.paridx

load('alldata.mat')
addpath(genpath('bads-master'));

nsubj = length(DATA);
nruns = 10;

nmodels = modelidx;
modelresults = cell(nsubj,nruns);
model = C_specifymodel(modelidx);

data = DATA{subjidx};
stim = STIM{subjidx};

% Fit parameters by subject, plug fitted parameters back into model 
pars_run = NaN(1,model.npars);
NLL_run = NaN;

model.init = unifrnd(model.plb, model.pub);  

[pars_run(1,:), NLL_run] = bads(@(par) C_modelpredictions(par, modelidx, stim, data), model.init, model.lb, model.ub, model.plb, model.pub);

[~, p_resp, false_d]       = C_modelpredictions(pars_run, modelidx, stim);

% Saving for this model and this subject
modelresults{subjidx,run}.par_est   = pars_run;
modelresults{subjidx,run}.NLL       = NLL_run;
modelresults{subjidx,run}.modelpred = p_resp;
modelresults{subjidx,run}.false_d   = false_d;

save(['results_' num2str(modelidx) '_sub_' num2str(subjidx) '_run_' num2str(run) '.mat'],'modelresults')

end
