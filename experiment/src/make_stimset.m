function stim = make_stimset(param)

nTrials     = param.nTrials;
ns          = param.ns;
ns          = ns';
sig_s       = param.sig_s;
sig         = param.sig;
R           = param.rad;
pC1         = param.pC1;

if rem(nTrials,length(ns)) > 0 && rem(nTrials,2) > 0
    fprintf('ERROR: set nTrials divisible by ns and 2')
end 

stim.feeder                       = zeros(nTrials,1);
stim.feeder(1:round(nTrials*pC1)) = 1;
stim.feeder                       = logical(stim.feeder);

stim.n                          = repmat(ns,[nTrials/length(ns),1]);
stim.n1(~stim.feeder,1)         = 0;
stim.n1(stim.feeder,1)          = binornd(stim.n(stim.feeder),0.5);
stim.n0                         = stim.n - stim.n1;

for i_trial = 1:nTrials
    goodDraw = 0;
    while goodDraw == 0
        [stim.x{i_trial}, stim.y{i_trial}, stim.sx{i_trial}, stim.sy{i_trial}, stim.x1{i_trial}, stim.y1{i_trial}] = step1_generative(stim.n0(i_trial), stim.n1(i_trial), 1, param);
        if sum(stim.x{i_trial} > param.rad) + sum(stim.y{i_trial} > param.rad) == 0
            goodDraw = 1;
        else
            goodDraw = 0;
        end
    end
end

%%%% randperm
iperm = randperm(nTrials);

stim.feeder = stim.feeder(iperm);
stim.n      = stim.n(iperm);
stim.n1     = stim.n1(iperm);
stim.n0     = stim.n0(iperm);
stim.x      = stim.x(iperm);
stim.y      = stim.y(iperm);
stim.sx     = stim.sx(iperm);
stim.sy     = stim.sy(iperm);
stim.x1     = stim.x1(iperm);
stim.y1     = stim.y1(iperm);
end