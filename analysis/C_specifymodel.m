function model = C_specifymodel(modelidx, paffmodel)

% paffmodel -- if paffmodel is 1, will use version of model with
% N-dependent p_aff

% "plus" denotes 4 criteria
% "star" denotes ?

% (In parentheses number of parameters)
% Model  52: (A+) 1 with 4 K's  --Bayes 
% Model  53: (A) 1 with false sigma_s and sigma -- false Bayes
% Model  54: (A) 1 with false p_affiliation -- false Bayes

% Model  2: (H) Global mean dist
% Model  3: (H) (Max) Local density-- Gaussian convolve
% Model  4: (H) Min distance
% Model  5: (H) Mean NN distance
% Model  14: (H) (Sum) Local density -- Gaussian convolve
% Model  16: (H) Fraction of points (z) below NNdist threshold (d)

% Model  39: (C+) 17 with 4 K's --Centroid add-one (BHC)-- only z-conditioned
% Model  40: (C+) 26 with 4 K's --Most likely z-- z estimated by maximizing marginal likelihood of z (marginal posterior is the same)

% Model  41: (D+) 32 with 4 K's --Max marginal likelihood mu-conditioned, z-optimal threshold conditioned
% Model  42: (D+) 33 with 4 K's --Max marginal posterior mu-conditioned, z-optimal threshold conditioned

% Model  43: (D+) 34 with 4 K's --Max local density -- mu-conditioned, z-optimal threshold conditioned
% Model  44: (D+) 35 with 4 K's --Mean density -- mu-conditioned, z-optimal threshold conditioned
% Model  45: (D+) 18 with 4 K's --Centroid add-one-- z-conditioned and mu-conditioned
% Model  46: (D+) 37 with 4 K's --Brute max joint likelihood
% Model  47: (D+) 38 with 4 K's --Brute max joint posterior

% Model  49: (B+) -- mu PME -- mu_conditioned, z-marginalized;
% Model  50: (B+) 29 with 4 K's --Max local density-- mu-conditioned, z-marginalized
% Model  51: (B+) 30 with 4 K's --Mean density -- mu-conditioned, z-marginalized




switch modelidx
      case {2, 3, 4, 5, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52}
        pars    = {'(1) k6',...
                   '(2) k9',...
                   '(3) k12',...
                   '(4) k15',...
                   '(5) sigma_d',...
                   '(6) lambda'};
        npars   = length(pars);
        lb      = [-100 -100 -100 -100  0   0 ];
        ub      = [ 100  100  100  100  100 1 ];
        plb     = [ -10 -10  -10  -10   0   0  ];
        pub     = [  10  10   10   10   10  0.1];
        
        if paffmodel
            pars{end+1} = '(7) p_aff6';
            pars{end+1} = '(8) p_aff9';
            pars{end+1} = '(9) p_aff12';
            pars{end+1} = '(10) p_aff15';
            
            npars   = length(pars);
            lb      = [lb 0 0 0 0];
            ub      = [ub 1 1 1 1];
            plb     = [plb 0 0 0 0];
            pub     = [pub 1 1 1 1];
        end
        
      case 53
        pars    = {'(1) k: decision criterion',...
                   '(2) sigma_d: decision noise',...
                   '(3) lambda: lapse rate',...
                   '(4) sigma_s: feeder spread around center',...
                   '(5) sigma: pigeon spread around feeder'};
        npars   = length(pars);
        lb      = [-inf, 0, 0, 0, 0];
        ub      = [inf, inf, 1, inf, inf];
        plb     = [-10, 0, 0, 0, 0];
        pub     = [10, 10, 0.2, 10, 10];
        
     case 54
        pars    = {'(1) k: decision criterion',...
                   '(2) sigma_d: decision noise',...
                   '(3) lambda: lapse rate',...
                   '(4) p_aff6',...
                   '(5) p_aff9',...
                   '(6) p_aff12',...
                   '(7) p_aff15'};
               
        npars   = length(pars);
        lb      = [-inf, 0, 0, 0, 0, 0, 0];
        ub      = [inf, inf, 1, 1, 1, 1, 1];
        plb     = [-10, 0, 0, 0, 0, 0, 0];
        pub     = [10, 10, 0.2, 1, 1, 1, 1];
        
        
     case 14
        pars    = {'(1) k6',...
                   '(2) k9',...
                   '(3) k12',...
                   '(4) k15',...
                   '(5) sigma_d',...
                   '(6) lambda'};
        npars   = length(pars);
        lb      = [0 0 0 0 0 0];
        ub      = [inf inf inf inf inf 1];
        plb     = [0 0 0 0 0 0];
        pub     = [3000 3000 3000 3000 10 0.2];
        
    case 6
        pars    = {'(1) sigma_d: decision noise',...
                   '(2) lambda: lapse rate'};
        npars   = length(pars);
        lb      = [0, 0];
        ub      = [inf, 1];
        plb     = [0, 0];
        pub     = [10, 1];
        
    case {15,23,24,25}
        pars    = {'(1) k: decision criterion',...
                   '(2) sigma_d: decision noise',...
                   '(3) lambda: lapse rate',...
                   '(4) like_thresh: likelihood threshold for z selection'};
        npars   = length(pars);
        lb      = [-inf, 0, 0, 0];
        ub      = [inf, inf, 1, 1];
        plb     = [-10, 0, 0, 0];
        pub     = [10, 10, 0.2, 0.2];
        
    case 16
        pars    = {'(1) d6',...
                   '(2) d9',...
                   '(3) d12',...
                   '(4) d15',...
                   '(5) z6',...
                   '(6) z9',...
                   '(7) z12',...
                   '(8) z15',...
                   '(9) lambda'};
        npars   = length(pars);
        lb      = [0, 0, 0, 0, 0, 0, 0, 0, 0];
        ub      = [inf, inf, inf, inf, 1, 1, 1, 1, 1];
        plb     = [0, 0, 0, 0, 0, 0, 0, 0, 0];
        pub     = [5, 5, 5, 5, 0.9, 0.9, 0.9, 0.9, 0.2];
        

end


switch modelidx
    case {2,3,4,5,14,16}
        cat     = 'H';
    case {1,53,54}
        cat     = 'A';
    case {27,28,29,30}
        cat     = 'B';
    case {17,26}
        cat     = 'C';
    case {18,32,33,34,35,37,38}
        cat     = 'D';
    case {6:13,19:22}
        cat     = 'X';
    case {31}
        cat     = 'Bstar';
    case {15,23:25,36}
        cat     = 'Dstar';
    case {39,40}
        cat     = 'Cplus';
    case {41,42,43,44,45,46,47}
        cat     = 'Dplus';
    case {48,49,50,51}
        cat     = 'Bplus';
    case {52}
        cat     = 'Aplus';
end

if paffmodel
    cat = [cat 'A'];
end

model.npars           = npars;
model.lb              = lb;
model.ub              = ub;
model.plb             = plb;
model.pub             = pub;
model.cat             = cat;
model.paff            = paffmodel;