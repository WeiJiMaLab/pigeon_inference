function generate_park_map(x0,y0,x1,y1,sx,sy,R,i_trial)

% Generates the park map of an example trial (to visualize multiple trials,
% loop)

% Input:
% x0,y0,x1,y1: matrix of coordinates of unaffiliated and affiliated pigeons
        % Size [nTrials, nPigeon0 or nPigeon1]
% sx,sy: matrix coordinates of the feeder
        % Size [nTrials, 1]
% R: park radius
% i_trial: the index of the trial 

% Output: opens figure plotting trial

    figd;
    scatter(x0(i_trial,:),y0(i_trial,:),100,'k','filled','LineWidth',3); hold on
    if(length(x1)>0)
        scatter(x1(i_trial,:),y1(i_trial,:),100,'b','filled','LineWidth',3); hold on
        scatter(sx(i_trial,:), sy(i_trial,:),100,'r','LineWidth',3); hold on
    end

    viscircles([0,0],R,'Color','k');

    title(['true (s_x,s_y) = ' num2str(sx(i_trial,:)) ', '  num2str(sy(i_trial,:))])
    axis square
    axis equal
    xlim([-10 10])
    ylim([-10 10])
    grid on
    xticks([0])
    yticks([0])
%     xlim([min(s) max(s)]);
%     ylim([min(s) max(s)]);
end