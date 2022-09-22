%% get responses
function [data, exitNow] = run_trial(param, stim, trialNum, data, isPract)

Screen('BlendFunction', param.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

exitNow = 0;
window = param.window;
windex = param.windex;

%% Draw the pigeons on the screen one-by-one
Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
data.timing.trialStart(trialNum) = Screen('Flip', window);
%WaitSecs(0.3)
%KbQueueFlush;

i_pigeon = length(stim.x{trialNum});
draw_pigeons(i_pigeon,stim.x{trialNum},-stim.y{trialNum},param);
Screen('DrawingFinished', param.window); 
Screen('Flip', param.window);
WaitSecs(0.4)


%%%%%%%%%%% start rating %%%%%%%%%%%

r1_text = ' ';

Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
DrawFormattedText(window, r1_text, 'center', 'center', param.fontColor, [], [], [], param.vSpacing);

WaitSecs(param.stimTime);
data.timing.r1_prompt(trialNum) = Screen('Flip',window);

% collect rating
[r1_key data.r1_RT(trialNum) data.timing.r1_time(trialNum)] = ... 
    recordValidKeys(data.timing.r1_prompt(trialNum), param.respDur_inSecs, param.keyboardNumber, param.r1_validKeys);

Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
Screen('Flip',window);

% assess response
id = 6;
rk = param.r1_validKeys;
switch r1_key
    case rk{1},         data.resp(trialNum) =  -1;   data.resp_feeder(trialNum) = 0; data.resp_conf(trialNum) = 1;
    case rk{2},         data.resp(trialNum) =  -2;   data.resp_feeder(trialNum) = 0; data.resp_conf(trialNum) = 2;
    case rk{3},         data.resp(trialNum) =  -3;   data.resp_feeder(trialNum) = 0; data.resp_conf(trialNum) = 3;
    case rk{4},         data.resp(trialNum) =  -4;   data.resp_feeder(trialNum) = 0; data.resp_conf(trialNum) = 4;
    case rk{5},         data.resp(trialNum) =  1;   data.resp_feeder(trialNum) = 1; data.resp_conf(trialNum) = 1;
    case rk{6},         data.resp(trialNum) =  2;   data.resp_feeder(trialNum) = 1; data.resp_conf(trialNum) = 2;
    case rk{7},         data.resp(trialNum) =  3;   data.resp_feeder(trialNum) = 1; data.resp_conf(trialNum) = 3;
    case rk{8},         data.resp(trialNum) =  4;   data.resp_feeder(trialNum) = 1; data.resp_conf(trialNum) = 4;
    case 'noanswer',    data.resp(trialNum) = -5;   data.resp_feeder(trialNum) = -5; data.resp_conf(trialNum) = 0;
    case 'invalid',     data.resp(trialNum) = -6;   data.resp_feeder(trialNum) = -6; data.resp_conf(trialNum) = 0;
    case 'cell',        data.resp(trialNum) = -7;   data.resp_feeder(trialNum) = -7; data.resp_conf(trialNum) = 0;
    case param.exitKey, data.resp(trialNum) = -8; exitNow = 1;
    otherwise,          data.resp(trialNum) = -9;   data.resp_feeder(trialNum) = -9; data.resp_conf(trialNum) = 0;
end

%% sort out response 

if exitNow, return; end

% sort out correctness
if data.resp_feeder(trialNum) == stim.feeder(trialNum)
    data.correct_feeder(trialNum) = 1;
    
    data.wager_val(trialNum) = data.resp_conf(trialNum);
    fbColor = param.green;
    wager_sign = '+';
    
% N/A if response < 1
elseif data.resp_feeder(trialNum) < 0
    data.correct_feeder(trialNum) = -1;

% otherwise, just incorrect
else
    data.correct_feeder(trialNum) = 0;
    
    data.wager_val(trialNum) = data.resp_conf(trialNum);
    fbColor = param.red;
    wager_sign = '-';
end

%% show feedback

switch data.resp_feeder(trialNum)
    case 0, resp_text = '"absent."';
    case 1, resp_text = '"present."';
    case -1, resp_text = '';
end

switch stim.feeder(trialNum)
    case 0, truth_text = 'The feeder was absent.';
    case 1, truth_text = 'The feeder was present.';
end

switch data.correct_feeder(trialNum)
    case 1,  FB_text1 = ['Correct! \n' wager_sign num2str(data.resp_conf(trialNum))]; FB_text2 = ['\n\nYou reported ' resp_text];
    case 0,  FB_text1 = ['Incorrect! \n' wager_sign num2str(data.resp_conf(trialNum))]; FB_text2 = ['\n\n' truth_text '\nYou reported ' resp_text];
    case -1, FB_text1 = ['No response was recorded.'];
end

if isPract
    DrawFormattedText(window, [FB_text1 FB_text2], param.screenXpixels*0.05, 'center', param.fontColor, [], [], [], param.vSpacing);
    DrawFormattedText(window, 'unaffiliated pigeons', param.screenXpixels*0.8, param.screenYpixels*0.4, param.white, [], [], [], param.vSpacing);
    DrawFormattedText(window, '(Press any key to continue)', 'center', param.screenYpixels*0.95, param.white, [], [], [], param.vSpacing);

    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    draw_pigeons(length(stim.x{trialNum}),stim.x{trialNum},-stim.y{trialNum},param); % Draw unaffiliated
    
    if stim.feeder(trialNum)
        DrawFormattedText(window, 'hidden feeder location', param.screenXpixels*0.8, param.screenYpixels*0.5, param.green, [], [], [], param.vSpacing);
        DrawFormattedText(window, 'affiliated pigeons', param.screenXpixels*0.8, param.screenYpixels*0.6, param.purple, [], [], [], param.vSpacing);
        draw_pigeons(length(stim.sx{trialNum}),stim.sx{trialNum},-stim.sy{trialNum},param, param.green); % Draw feeder
        draw_pigeons(length(stim.x1{trialNum}),stim.x1{trialNum},-stim.y1{trialNum},param, param.purple); % Draw affiliated
    end
    
    Screen('DrawingFinished', param.window); 
    Screen('Flip', param.window);
    
    KbStrokeWait
        
    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    Screen('Flip', window);
    WaitSecs(0.5);
else
    DrawFormattedText(window, FB_text1, 'center', param.screenYpixels*0.4, fbColor, [], [], [], param.vSpacing);
    DrawFormattedText(window, FB_text2, 'center', 'center', param.fontColor, [], [], [], param.vSpacing);
    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    Screen('Flip', window);
    WaitSecs(0.4);
end

Screen('BlendFunction', param.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
data.timing.FBoffset(trialNum) = Screen('Flip', window);
data.timing.trialEnd(trialNum) = GetSecs;
data.timing.trialDur(trialNum) = data.timing.trialEnd(trialNum) - data.timing.trialStart(trialNum);
