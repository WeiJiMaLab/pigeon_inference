%% draw text screen
function [goBack] = draw_demo_text(i_screen,param,stim)
window = param.window;
commandwindow; %makes it so characters typed do show up in the command window

leftKey = KbName('LeftArrow');
escKey = KbName('ESCAPE');
spaceKey = KbName('space');
i_trial = 1;

Screen('FillRect',window, param.white);
trials_C1 = find(stim.feeder);
trials_C0 = find(~stim.feeder);

if param.pC1 == 0.5
    prior_text = 'HALF (50%)';
elseif param.pC1 == 0.333
    prior_text = 'ONE THIRD (33%)';
end

switch i_screen
    case 1
    %Introduce circular park
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Welcome to Paranormal Pigeon Park!\nYou''ll notice that the park boundaries are a circle.','center','center', param.white);
    
    case 2
    %Introduce center cross
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    Screen('BlendFunction', param.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    DrawFormattedText(param.window,'The center of the park is always given by this grey cross. ','center', param.screenYpixels*0.45, param.white);
    
    case 3
    %Introduce pigeons
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    draw_pigeons(length(stim.x{trials_C1(1)}),stim.x{trials_C1(1)},-stim.y{trials_C1(1)},param);
    DrawFormattedText(param.window,'Say hello to the pigeons! \nEach pigeon is represented by a white dot.','center', param.screenYpixels*0.45, param.blue);
    Screen('DrawingFinished', param.window);
    
    case 4
    %Introduce C=0
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'On HALF (50%) of the days, there is NO PIGEON FEEDER. \n\nThis means pigeons are scattered randomly within the circle','center', 'center', param.white);
    
    case 5
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'With NO FEEDER, the pigeons are all \nUNAFFILIATED PIGEONS. \n\n\nThese pigeons pick their location completely at random.','center', 'center', param.blue);
    
    case 6
    %Example C=0
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Here are some examples of what the park looks like \nwith NO PIGEON FEEDER.','center', 'center', param.white);
    
    case 7
    %C=0 examples
    while i_trial <= length(trials_C0) && i_trial > 0
        trialNum = trials_C0(i_trial);
        my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
        draw_pigeons(length(stim.x{trialNum}),stim.x{trialNum},-stim.y{trialNum},param);
        DrawFormattedText(param.window,'<-- Back                                                             Next -->','center', param.screenYpixels*0.95, param.white);
        Screen('DrawingFinished', param.window);
        Screen('Flip', param.window); 
        [~, keyCode] = KbStrokeWait;
        if find(keyCode) == leftKey && i_screen > 1
            i_trial = i_trial -1;
        elseif find(keyCode) == escKey
            sca
        else
            i_trial = i_trial + 1;
        end
    end
    if i_trial > length(trials_C0)
        goBack = 0;
        return
    elseif i_trial == 0
        goBack = 1;
        return
    end
    
    case 8
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,['Say hello to the PIGEON FEEDER! \n\n\nShe''s a lonely old GHOST who appears in the \npark on ' prior_text ' of the days.'],'center', 'center', param.red);
    draw_pigeons(length(stim.sx{1}),stim.sx{1},-stim.sy{1},param, param.green); % Draw feeder

    case 9
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Legend has it that she died here long ago. \n\nShe typically appears NEAR THE CENTER of the park, \nbut she tends to stray.','center', 'center', param.white);
    
    case 10
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Here are some examples of where the ghostly \nPIGEON FEEDER has been spotted in the past!','center', 'center', param.red);
    
    case 11
    %Pigeon feeder location examples
    while i_trial <= length(trials_C1) && i_trial > 0
        trialNum = trials_C1(i_trial);
        
        DrawFormattedText(window, 'pigeon feeder sighting', param.screenXpixels*0.8, param.screenYpixels*0.5, param.green, [], [], [], param.vSpacing);
        draw_pigeons(length(stim.sx{trialNum}),stim.sx{trialNum},-stim.sy{trialNum},param, param.green); % Draw feeder
        DrawFormattedText(param.window,'<-- Back                                                             Next -->','center', param.screenYpixels*0.95, param.white);
        Screen('DrawingFinished', param.window);
        Screen('Flip', param.window);
        
        [~, keyCode] = KbStrokeWait;
        if find(keyCode) == leftKey && i_screen > 1
            i_trial = i_trial-1;
        elseif find(keyCode) == spaceKey
            i_trial = i_trial;
        elseif find(keyCode) == escKey
            sca
        else
            i_trial = i_trial + 1;
        end
    end
    
    if i_trial > length(trials_C1)
        goBack = 0;
        return
    elseif i_trial == 0
        goBack = 1;
        return
    end
    
    
    case 12
    %Explain C=1
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,['On exactly ' prior_text ' of all days, she visits the park.'],'center', 'center', param.red);
    
    case 13
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'When the ghostly PIGEON FEEDER is present, \npigeons tend to cluster around her location.','center', 'center', param.white);

    case 14
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'... But even pigeons can''t always sense her presence!','center', 'center', param.white);

    case 15
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Even when the FEEDER is PRESENT, there''s only a \n50% chance that a given pigeon will be AFFILIATED with her.','center', 'center', param.red);
    
    case 16
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Here are some examples of what the park looks like \nwhen the FEEDER is PRESENT.','center', 'center', param.white);

    case 17
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'AFFILIATED pigeons tend to cluster near the pigeon feeder. \n\nPay attention to the SPREAD of AFFILIATED PIGEONS \naround the feeder!','center', 'center', param.red);

    case 18
    %C=1 examples
    while i_trial <= length(trials_C1) && i_trial > 0
        trialNum = trials_C1(i_trial);
        my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
        draw_pigeons(length(stim.x{trialNum}),stim.x{trialNum},-stim.y{trialNum},param);
        DrawFormattedText(param.window,'<-- Back          (TOGGLE THE SPACE BAR TO REVEAL THE FEEDER)          Next -->','center', param.screenYpixels*0.95, param.white);
        Screen('Flip', param.window);
        [~, keyCode] = KbStrokeWait;
        if find(keyCode) == leftKey && i_screen > 1
            i_trial = i_trial-1; continue
        elseif find(keyCode) == escKey
            sca
        end
        
        DrawFormattedText(window, 'unaffiliated pigeons', param.screenXpixels*0.8, param.screenYpixels*0.4, param.white, [], [], [], param.vSpacing);
        DrawFormattedText(window, 'hidden feeder location', param.screenXpixels*0.8, param.screenYpixels*0.5, param.green, [], [], [], param.vSpacing);
        DrawFormattedText(window, 'affiliated pigeons', param.screenXpixels*0.8, param.screenYpixels*0.6, param.purple, [], [], [], param.vSpacing);
        DrawFormattedText(param.window,'<-- Back          (TOGGLE THE SPACE BAR TO REVEAL THE FEEDER)          Next -->','center', param.screenYpixels*0.95, param.white);
        draw_pigeons(length(stim.x{trialNum}),stim.x{trialNum},-stim.y{trialNum},param); % Draw unaffiliated
        draw_pigeons(length(stim.sx{trialNum}),stim.sx{trialNum},-stim.sy{trialNum},param, param.green); % Draw feeder
        draw_pigeons(length(stim.x1{trialNum}),stim.x1{trialNum},-stim.y1{trialNum},param, param.purple); % Draw affiliated
        Screen('DrawingFinished', param.window);
        Screen('Flip', param.window);
        
        [~, keyCode] = KbStrokeWait;
        if find(keyCode) == leftKey && i_screen > 1
            i_trial = i_trial-1;
        elseif find(keyCode) == spaceKey
            i_trial = i_trial;
        elseif find(keyCode) == escKey
            sca
        else
            i_trial = i_trial + 1;
        end
    end
    
    if i_trial > length(trials_C1)
        goBack = 0;
        return
    elseif i_trial == 0
        goBack = 1;
        return
    end
    
    case 19
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Your mission at paranormal pigeon park is to tell apart the days \non which the ghostly PIGEON FEEDER is \n\nABSENT (left-hand) vs PRESENT (right-hand).','center', 'center', param.white);

    case 20
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Please use the full range of confidence values when you respond!','center', 'center', param.white);
    
    case 21
    my_circle(window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    DrawFormattedText(param.window,'Time for some practice rounds!','center', 'center', param.white);
end

DrawFormattedText(param.window,'<-- Back                                                             Next -->','center', param.screenYpixels*0.95, param.white);
Screen('Flip', param.window);
[~, keyCode] = KbStrokeWait;

if find(keyCode) == leftKey && i_screen > 1
    goBack = 1;
elseif find(keyCode) == escKey
    sca
else
    goBack = 0;
end

end
    