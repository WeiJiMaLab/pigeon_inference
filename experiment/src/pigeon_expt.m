clear
clc
close all

%%
%%%% BEGIN PARAMETERS %%%%
isPract = 1;

block_num = '10';
subjectID = 'JL';
nTrials = 200;
ns = [6 9 12 15];
pC1 = 0.5;

if isPract
    block_num = 'pract';
    nTrials   = 20;
    ns        = [15];
    pC1       = 0.5;
end

sig_s = 2;
sig = 2;

%%%% END PARAMETERS %%%%

%% Initialize experiment
addpath(genpath('/Users/jennlauralee/Google Drive/WedMock/Causal Inference/Pigeon_expt1'))
Screen('Preference', 'SkipSyncTests', 1);
HideCursor();

% PTB default setup
PsychDefaultSetup(2);
commandwindow; %makes it so characters typed do show up in the command window

%params = set_params(block_num, sub_ID, nTrials, ns, sig_s, sig, pC1);
param = set_params(block_num, subjectID, nTrials, ns, sig_s, sig, pC1);

% Make or load stimulus x and y
stim = make_stimset(param);

%% Draw instructions
draw_instructions(param);

%% Run through trials
data=[];
tic
for trialNum = 1:param.nTrials
    [data, exitNow] = run_trial(param, stim, trialNum, data, isPract);
    if exitNow
        save([param.sub_ID '_iiiblock' param.block_num '_SC_' datestr(now, 'yyyymmdd-HHMMSS') '.mat'], 'stim', 'data', 'param');
        ShowCursor();
        sca;
    end
end
toc

%% save variables
save([param.sub_ID '_block' param.block_num '_SC_' datestr(now, 'yyyymmdd-HHMMSS') '.mat'], 'stim', 'data', 'param');

%% Done block
Screen('FillRect',param.window, param.black);
Screen('Flip', param.window);
DrawFormattedText(param.window,'Nice work! Break time! \n\nPlease find Jenn.','center','center', param.white);
Screen('Flip', param.window);
KbStrokeWait;

ShowCursor();
sca;