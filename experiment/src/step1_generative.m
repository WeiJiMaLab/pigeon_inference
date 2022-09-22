%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Generative model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y, sx, sy, x1, y1] = step1_generative(n0, n1, ntrials, param)
% Generates x and y pair coordinates of pigeon locations

% Input:
% n0: number of unaffiliated pigeons
% n1: number of affiliated pigeons
% ntrials: number of simulated trials
% param.sig_s: how far from center does the feeder spawn?
% param.sig: how far from feeder do pigeons spawn?
% param.rad: radius of the circle

% Output: x and y pair coordinate of pigeon locations (each ntrials x 1)
%% set parameters
if nargin > 3
    sig_s       = param.sig_s; % how far from center does the feeder spawn?
    sig         = param.sig; %how far from feeder do pigeons spawn?
    R           = param.rad; %radius of the circle
else
    sig_s       = 2;
    sig         = 2;
    R         = 10;
end

mu_s        = 0; % center of the park around which feeders congregate
s           = -R:0.1:R;

%% generate feeder locations (s) -- in 2AFC case there is only 1 feeder
sx = normrnd(mu_s, sig_s, [ntrials,1]); %generate true x-coordinates of the feeders; K-1 because 1 is null
sy = normrnd(mu_s, sig_s, [ntrials,1]); %generate true y-coordinates of the feeders

%% Draw from uniform distribution (limited by the circular radius of the screen)
theta_draw  = 2*pi.*rand(ntrials,n0);
r_draw      = R.*sqrt(rand(ntrials,n0));
x0          = (r_draw).*cos(theta_draw);
y0          = (r_draw).*sin(theta_draw);

%% Draw from causal distribution
x1 = [];
y1 = [];
if n1 > 0
    for ii = 1:n1
        x1 = [x1 normrnd(sx, sig)]; %x-coordinates of each pigeon
        y1 = [y1 normrnd(sy, sig)]; %y-coordinates of each pigeon
    end
end

x = [x0 x1];
y = [y0 y1];

if n0 == 0 & n1 == 0
    x = nan(ntrials,1);
    y = nan(ntrials,1);
end

%generate_park_map(x0,y0,x1,y1,sx,sy,R,1)

end