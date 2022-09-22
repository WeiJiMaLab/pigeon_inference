%% C_modelfitting

function C_modelfitting_cluster(i_job)
% make each batch its own model

addpath(genpath(pwd))

subvec = [1:10];
runvec = [1:100];
paffvec = [1];
modelidx = [9];

designmat = combvec(subvec,runvec,paffvec,modelidx)';

subjidx   = designmat(i_job,1);
run       = designmat(i_job,2);
paffmodel = designmat(i_job,3);
modelidx  = designmat(i_job,4);

% Todo: save general settings above subject level in superstructure e.g.
% modelresults{modelidx}. and modelresults{modelidx}.bysubject{subjidx}.paridx

load('alldata.mat')
addpath(genpath('bads-master'));
addpath(genpath('helper-functions'));

nsubj = length(DATA);
nruns = 100;

nmodels = modelidx;
modelresults = cell(nsubj,nruns);
model = C_specifymodel_v2(modelidx, paffmodel);

data = DATA{subjidx};
stim = STIM{subjidx};

switch modelidx
    case {9,10} %{46,47}
        load(['precomputed_LLR_forallz_S' num2str(subjidx) '.mat']);
        precomputed_LLR.like = precomputed_LLR_forallz.like{subjidx};
        precomputed_LLR.post = precomputed_LLR_forallz.post{subjidx};
end

% Fit parameters by subject, plug fitted parameters back into model 
pars_run = NaN(1,model.npars);
NLL_run = NaN;

model.init = unifrnd(model.plb, model.pub);  

[pars_run(1,:), NLL_run] = bads(@(par) C_modelpredictions_v2(par, modelidx, stim, data, paffmodel), model.init, model.lb, model.ub, model.plb, model.pub);

%[~, p_resp, false_d]       = C_modelpredictions(pars_run, modelidx, stim, [], paffmodel, precomputed_LLR);
[~, p_resp]       = C_modelpredictions_v2(pars_run, modelidx, stim, [], paffmodel);


% Saving for this model and this subject
modelresults{subjidx,run}.par_est   = pars_run;
modelresults{subjidx,run}.NLL       = NLL_run;
modelresults{subjidx,run}.modelpred = p_resp;
%modelresults{subjidx,run}.false_d   = false_d;

if paffmodel
    pafftext = '_paff';
else
    pafftext = '';
end

save(['results_' num2str(modelidx) pafftext '_sub_' num2str(subjidx) '_run_' num2str(run) '.mat'],'modelresults')

end
