function viz_park_map(x0,y0,x1,y1,R,i_trial)

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

if length(x0>0)
    scatter(x0(i_trial,:),y0(i_trial,:),'k','filled','LineWidth',3); hold on
end
if length(x1>0)
    scatter(x1(i_trial,:),y1(i_trial,:),'b','filled','LineWidth',3); hold on
end

%scatter(sx(i_trial,:), sy(i_trial,:),'r','LineWidth',3); hold on

viscircles([0,0],R,'Color','k');

%title(['true (s_x,s_y) = ' num2str(sx(i_trial,:)) ', '  num2str(sy(i_trial,:))])
axis square
axis equal
grid on

end