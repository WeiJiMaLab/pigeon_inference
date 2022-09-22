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
%KbQueueFlush;

i_pigeon = length(stim.x{trialNum});
draw_pigeons(i_pigeon,stim.x{trialNum},-stim.y{trialNum},param);
Screen('DrawingFinished', param.window); 
Screen('Flip', param.window);
WaitSecs(0.4)

%%%%%%%%%%% start rating %%%%%%%%%%%

r1_text = 'Was there a feeder?\n\n F = no        J = yes';

DrawFormattedText(window, r1_text, 'center', param.screenYpixels*0.4, param.fontColor, [], [], [], param.vSpacing);
Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim

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
    case rk{1},         data.resp_feeder(trialNum) =  0;  id = 0; resp_feeder = 0;
    case rk{2},         data.resp_feeder(trialNum) =  1;  id = 1; resp_feeder = 1;
    case 'noanswer',    data.resp_feeder(trialNum) = -1;  id = -1; resp_feeder = -1;
    case 'invalid',     data.resp_feeder(trialNum) = -2;  id = -2; resp_feeder = -1;
    case 'cell',        data.resp_feeder(trialNum) = -3;  id = -3; resp_feeder = -1;
    case param.exitKey, data.resp_feeder(trialNum) = -4; exitNow = 1;
    otherwise,          data.resp_feeder(trialNum) = -5;  id = -5; resp_feeder = -1;
end

% give 'too slow' warning if necessary
if strcmp(r1_key,'noanswer')
    data.timing.r1_tooSlow(trialNum) = 1;
    DrawFormattedText(window, 'Too slow!', 'center', param.screenYpixels*0.4, param.fontColor)
    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    Screen('Flip', window);
    WaitSecs(2);
end

%%%%%%%%%%% confidence query %%%%%%%%%%%

r2_text = 'Confidence?\n\n 1 key --> low confidence \n4 key --> high confidence';
DrawFormattedText(window, r2_text, 'center', param.screenYpixels*0.4, param.fontColor, [], [], [], param.vSpacing);
Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim

WaitSecs(.4);
data.timing.r2_prompt(trialNum) = Screen('Flip',window);

% collect beta estimation
if data.resp_feeder(trialNum) == 0
    r2_validKeys = param.r2_validKeys0;
elseif data.resp_feeder(trialNum) == 1
    r2_validKeys = param.r2_validKeys1;
else
    r2_validKeys = param.r2_validKeys;
end

rk = r2_validKeys;
[r2_key data.r2_RT(trialNum) data.timing.r2_time(trialNum)] = ... 
    recordValidKeys(data.timing.r2_prompt(trialNum), param.respDur_inSecs, param.keyboardNumber, r2_validKeys);

% assess response
id = 6;

switch r2_key
    case rk{1},         data.resp_conf(trialNum) =  1;   id = 4;  
    case rk{2},         data.resp_conf(trialNum) =  2;   id = 3; 
    case rk{3},         data.resp_conf(trialNum) =  3;   id = 2;   
    case rk{4},         data.resp_conf(trialNum) =  4;   id = 1;
    case 'noanswer',    data.resp_conf(trialNum) = -1;
    case 'invalid',     data.resp_conf(trialNum) = -2;
    case 'cell',        data.resp_conf(trialNum) = -3;
    case param.exitKey, data.resp_conf(trialNum) = -4; exitNow = 1;
    otherwise,          data.resp_conf(trialNum) = -5;
end

% give 'too slow' warning if necessary
if strcmp(r2_key,'noanswer')
    data.timing.r2_tooSlow(trialNum) = 1;    
    DrawFormattedText(window, 'Too slow!', 'center', 'center', param.fontColor)
    Screen('Flip', window);
    WaitSecs(2);
end


%% sort out response 

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
    
    data.wager_val(trialNum) = - data.resp_conf(trialNum);
    fbColor = param.red;
    wager_sign = '-';
end

if exitNow, return; end


%% show feedback

switch resp_feeder
    case 0, resp_text = '"absent."';
    case 1, resp_text = '"present."';
    case -1, resp_text = '';
end

switch stim.feeder(trialNum)
    case 0, truth_text = 'The feeder was absent.';
    case 1, truth_text = 'The feeder was present.';
end

switch data.correct_feeder(trialNum)
    case 1,  FB_text = ['Correct! ' wager_sign num2str(data.resp_conf(trialNum)) '\n\nYou reported ' resp_text];
    case 0,  FB_text = ['Incorrect! ' wager_sign num2str(data.resp_conf(trialNum)) '\n\n' truth_text '\nYou reported ' resp_text];
    case -1, FB_text = ['No response was recorded.'];
end

if isPract
    DrawFormattedText(window, FB_text, param.screenXpixels*0.05, 'center', fbColor, [], [], [], param.vSpacing);
    DrawFormattedText(window, 'unaffiliated pigeons', param.screenXpixels*0.8, param.screenYpixels*0.4, param.white, [], [], [], param.vSpacing);
    DrawFormattedText(window, '(Press any key to continue)', 'center', param.screenYpixels*0.95, param.white, [], [], [], param.vSpacing);

    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    draw_pigeons(length(stim.x{trialNum}),stim.x{trialNum},-stim.y{trialNum},param); % Draw unaffiliated
    
    if stim.feeder(trialNum)
        DrawFormattedText(window, 'hidden feeder location', param.screenXpixels*0.8, param.screenYpixels*0.5, param.yellow, [], [], [], param.vSpacing);
        DrawFormattedText(window, 'affiliated pigeons', param.screenXpixels*0.8, param.screenYpixels*0.6, param.purple, [], [], [], param.vSpacing);
        draw_pigeons(length(stim.sx{trialNum}),stim.sx{trialNum},-stim.sy{trialNum},param, param.yellow); % Draw feeder
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
    DrawFormattedText(window, FB_text, 'center', param.screenYpixels*0.4, param.fontColor, [], [], [], param.vSpacing);
    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    Screen('Flip', window);
    WaitSecs(0.3);
end

Screen('BlendFunction', param.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
data.timing.FBoffset(trialNum) = Screen('Flip', window);
data.timing.trialEnd(trialNum) = GetSecs;
data.timing.trialDur(trialNum) = data.timing.trialEnd(trialNum) - data.timing.trialStart(trialNum);
%ShowCursor(); %shows the cursor
