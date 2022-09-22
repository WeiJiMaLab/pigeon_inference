%% demo

function practice(param)

%%
pract=[];
for trialNum = 1:param.nTrials
    [pract, exitNow] = run_trial(param, stim, trialNum, pract);
    if exitNow
        save([sub_ID '_iiiblock' block_num '_SC_' datestr(now, 'yyyymmdd-HHMMSS') '.mat'], 'stim', 'data', 'param');
        sca;
    end
end