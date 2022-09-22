function draw_pigeons(i_pigeon,mux,muy,param,color)
    
    if nargin < 5
        color = param.white;
    end

    dotsize = 9; % dot diameter (pixels) 
    [param.midW, param.midH] = getScreenMidpoint(param.windex); %midpoints

    x = mux(1:i_pigeon).*param.scale; 
    y = muy(1:i_pigeon).*param.scale;

    Screen('BlendFunction', param.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    Screen('DrawLines', param.window, param.allCoords, 2, param.grey, [param.midW param.midH], 2);
    my_circle(param.window,param.grey, param.midW, param.midH, param.rad*param.scale, param.grey,1); % Draw rim
    
    Screen('DrawDots', param.window, [x ; y],dotsize,color,[param.midW param.midH],1); % Draw dots
end