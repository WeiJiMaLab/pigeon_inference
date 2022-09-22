%% C_modelfitting
clear
close all

paffmodel = 0;
for modelidx = [1:2,4:13]
    % Todo: save general settings above subject level in superstructure e.g.
    % modelresults{modelidx}. and modelresults{modelidx}.bysubject{subjidx}.paridx

    load('alldata.mat')

    nsubj = length(DATA);
    nruns = 10;

    nmodels = modelidx;
    modelresults = cell(nmodels, nsubj);
    modelidx
    model = C_specifymodel_v2(modelidx,paffmodel);
    
    for subjidx = 1:nsubj
        tic
        subjidx
        data = DATA{subjidx};
        stim = STIM{subjidx};
        
        % Fit parameters by subject, plug fitted parameters back into model 
        pars_run = NaN(nruns,model.npars);
        NLL_run = NaN(1,nruns);
             
        for run = 1:nruns
            run
            model.init = unifrnd(model.plb, model.pub);
            %model.init = [-log(2^6)+7.8029     -log(2^9)+8.2296    -log(2^12)+8.5191    -log(2^15)+8.5311    1.8769         0.0020 0.5 0.5 0.5 0.5];
            
            %C_modelpredictions(model.init, modelidx, stim, data, 1)
            OPTIONS = bads('defaults');
            [pars_run(run,:), NLL_run(run)] = bads(@(par) C_modelpredictions_v2(par, modelidx, stim, data, paffmodel), model.init, model.lb, model.ub, model.plb, model.pub, [], OPTIONS);
        end
        NLL_run
        [NLL, runidx] = min(NLL_run);
        maxdiffNLL = max(NLL_run) - min(NLL_run);
        par_est = pars_run(runidx,:);
        [~, p_resp]       = C_modelpredictions_v2(par_est, modelidx, stim, [], paffmodel);
        
        % Saving for this model and this subject
        modelresults{modelidx,subjidx}.par_est   = par_est;
        modelresults{modelidx,subjidx}.NLL       = NLL;
        modelresults{modelidx,subjidx}.maxdiffNLL= maxdiffNLL;
        modelresults{modelidx,subjidx}.NLL_run   = NLL_run;
        modelresults{modelidx,subjidx}.modelpred = p_resp;
        
        toc
    end

%% save in general superstructure
for subjidx = 1:nsubj
    Modelresults{modelidx}.par_est(subjidx,:)   = modelresults{modelidx,subjidx}.par_est;
    Modelresults{modelidx}.NLL(subjidx,:)       = modelresults{modelidx,subjidx}.NLL;
    Modelresults{modelidx}.maxdiffNLL(subjidx)  = modelresults{modelidx,subjidx}.maxdiffNLL;
    Modelresults{modelidx}.NLL_run(subjidx,:)   = modelresults{modelidx,subjidx}.NLL_run;
    Modelresults{modelidx}.modelpred(subjidx,:) = modelresults{modelidx,subjidx}.modelpred;
end

if paffmodel
    pafftext = ['_paff'];
else
    pafftext = [''];
end
resultspath = '/Users/jennlauralee/GitHub Repos/pigeons/Analysis/results/';
save([resultspath 'results_' num2str(modelidx) pafftext '.mat'],'DATA','STIM','Modelresults')

%clear DATA STIM Modelresults modelresults
end

%D_plotdatawithfits