%% demo video

clear
clc
close all

% Set prior value
pC1 = 0.5;

addpath([pwd '/functions']);
Screen('Preference', 'SkipSyncTests', 1);
HideCursor();

% PTB default setup
%PsychDefaultSetup(2);
commandwindow; %makes it so characters typed do show up in the command window

%params = set_params(block_num, sub_ID, nTrials, ns, sig_s, sig);
param = set_params('test', 'testsub', 20, [15], 2, 2, pC1);
window = param.window;

% Make or load stimulus x and y
%load(['/Users/jennlauralee/Google Drive/WedMock/Causal Inference/stimulus examples/5n0_5n1.mat'])
stim = make_stimset(param);

lastScreen = 21;
i_screen = 1;
while i_screen <= lastScreen
    goBack = draw_demo_text(i_screen,param,stim);
    if goBack
        i_screen = i_screen - 1;
    else
        i_screen = i_screen + 1;  
    end
end

ShowCursor();
sca