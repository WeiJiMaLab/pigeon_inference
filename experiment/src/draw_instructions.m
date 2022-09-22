function draw_instructions(param)

Screen('FillRect',param.window, param.black);
Screen('Flip', param.window);
DrawFormattedText(param.window,'Pigeons will appear on the screen. You will be asked to report whether or not there is likely a pigeon feeder.\n\n\n\n (Press any key to continue.)','center','center', param.white);
Screen('Flip', param.window);
KbStrokeWait;

Screen('FillRect',param.window, param.black);
Screen('Flip', param.window);
DrawFormattedText(param.window,['When you respond, you will report a particular CONFIDENCE LEVEL. \n\n\n\n '...
    'If you respond Feeder ABSENT, use your left hand keys [F D S A] to report your confidence on a scale of 1-4. \n\n' ...
    'If you respond Feeder PRESENT, use your right hand keys [J K L ;] to report your confidence on a scale of 1-4. \n\n\n\n (Press any key to continue.)'],'center','center', param.white);
Screen('Flip', param.window);
KbStrokeWait;

Screen('FillRect',param.window, param.black);
Screen('Flip', param.window);
DrawFormattedText(param.window,[ 'For both hands: Pinky key = HIGH confidence. Index key = LOW confidence.' ... 
    '\n\n\n\n (Press any key to continue.)'],'center','center', param.white);
Screen('Flip', param.window);
KbStrokeWait;

end