function param = set_params(block_num, sub_ID, nTrials, ns, sig_s, sig, pC1)

if nargin == 0
    % Define default (practice) trial parameters
    param.nTrials       = 10;
    param.block_num     = 'test';
    param.sub_ID        = 'test';
    param.ns            = [15];
    param.sig_s         = 2;
    param.sig           = 2;
    param.pC1           = 0.5;
else
    param.nTrials       = nTrials;
    param.block_num     = block_num;
    param.sub_ID        = sub_ID;
    param.ns            = ns;
    param.sig_s         = sig_s;
    param.sig           = sig;
    param.pC1           = pC1;
end

param.rad     = 10;
param.scale   = 35;

% Get screen numbers, use max screen number
screens = Screen('Screens');
param.windex = max(screens);

% Define black and white
param.white = WhiteIndex(param.windex);
param.black = BlackIndex(param.windex);
param.grey = [0.2 0.2 0.2];
param.red = [1 0 0];
param.green = [0 1 0];
param.blue = [0 1 1];
param.purple =  [1 0 1];
param.yellow = [1 1 0];

% Define fixation cross
fixCrossDimPix = 10;
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
param.allCoords = [xCoords; yCoords];

% Open black screen
[param.window, windowRect] = PsychImaging('OpenWindow', param.windex, param.black);

% Define visual parameters
% text
param.BGcolour = param.grey;
param.font = 'Helvetica';
param.fontColor = 255;
param.textSize = 30;
param.sx = 200;
param.vSpacing = 1.4;
param.stimTime = 0.4;

% Define response parameters
param.keyboardNumber = getKeyboardNumber;
param.respDur_inSecs = 100; % time allowed for each response before triggering next trial automatically
param.exitKey        = 'ESCAPE';
param.readyKey       = 'Return';
param.r1_validKeys   = {'f' 'd' 's' 'a' 'j' 'k' 'l' ';:' param.exitKey};

[param.midW, param.midH] = getScreenMidpoint(param.windex);
[param.screenXpixels, param.screenYpixels] = Screen('WindowSize', param.window); % Get # pixels
