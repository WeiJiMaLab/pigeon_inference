function model = C_specifymodel_v2(modelidx, paffmodel)

switch modelidx
    
     case 1
        pars    = {'(1) k: decision criterion',...
                   '(2) sigma_d: decision noise',...
                   '(3) lambda: lapse rate'};
               
        npars   = length(pars);
        lb      = [-inf, log(0), 0];
        ub      = [inf, log(inf), 1];
        plb     = [-10, log(0.01), 0];
        pub     = [10, log(10), 0.2];
        
         if paffmodel
            pars{end+1} = '(4) p_aff6';
            pars{end+1} = '(5) p_aff9';
            pars{end+1} = '(6) p_aff12';
            pars{end+1} = '(7) p_aff15';
            
            npars   = length(pars);
            lb      = [lb 0 0 0 0];
            ub      = [ub 1 1 1 1];
            plb     = [plb 0 0 0 0];
            pub     = [pub 1 1 1 1];
         end
        
     case {2,4,5,6,7,8,9,10,11,12,13}
        pars    = {'(1) k6',...
                   '(2) k9',...
                   '(3) k12',...
                   '(4) k15',...
                   '(5) sigma_d',...
                   '(6) lambda'};
        npars   = length(pars);
        lb      = [-inf -inf -inf -inf  log(0)   0 ];
        ub      = [ inf  inf  inf  inf  log(inf) 1 ];
        plb     = [ -10 -10  -10  -10   log(0.01)   0  ];
        pub     = [  10  10   10   10   log(10)  0.1];
        
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
        
    case 3
        pars    = {'(1) k: decision criterion',...
                   '(2) sigma_d: decision noise',...
                   '(3) lambda: lapse rate',...
                   '(4) sigma_s: feeder spread around center',...
                   '(5) sigma: pigeon spread around feeder'};
        npars   = length(pars);
        lb      = [-inf, log(0), 0, log(0), log(0)];
        ub      = [inf, log(inf), 1, log(inf), log(inf)];
        plb     = [-10, log(0.01), 0, log(0.01), log(0.01)];
        pub     = [10, log(10), 0.2, log(10), log(10)];
        
    case 14
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
        lb      = [log(0), log(0), log(0), log(0), 0, 0, 0, 0, 0];
        ub      = [log(inf), log(inf), log(inf), log(inf), 1, 1, 1, 1, 1];
        plb     = [log(0.01), log(0.01), log(0.01), log(0.01), 0, 0, 0, 0, 0];
        pub     = [log(5), log(5), log(5), log(5), 0.9, 0.9, 0.9, 0.9, 0.2];
        
     case {15,16,17,18,19}
        pars    = {'(1) k6',...
                   '(2) k9',...
                   '(3) k12',...
                   '(4) k15',...
                   '(5) sigma_d',...
                   '(6) lambda'};
        npars   = length(pars);
        lb      = [log(0) log(0) log(0) log(0) log(0) 0];
        ub      = [log(inf) log(inf) log(inf) log(inf) log(inf) 1];
        plb     = [log(0.01) log(0.01) log(0.01) log(0.01) log(0.01) 0];
        pub     = [log(3000) log(3000) log(3000) log(3000) log(10) 0.2];
end


switch modelidx
    case {1,2,3}
        cat     = 'A';
    case {4,5,6}
        cat     = 'B';
    case {7,8}
        cat     = 'C';
    case {9,10,11,12,13}
        cat     = 'D';
    case {14,15,16,17,18,19}
        cat     = 'H';
end

if paffmodel
    cat = [cat 'P'];
end

model.npars           = npars;
model.lb              = lb;
model.ub              = ub;
model.plb             = plb;
model.pub             = pub;
model.cat             = cat;
model.paff            = paffmodel;