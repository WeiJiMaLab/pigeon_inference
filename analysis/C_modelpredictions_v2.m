function [NLL, p_RF, false_d] = C_modelpredictions_v2(pars, modelidx, stim, data, paffmodel, precomputed_LLR)

% Fitting mode or model prediction mode?
nodata      = isempty(data);   

switch modelidx
    case 1 % (A) bayes 
        k       = pars(1);
        sigma_d = exp(pars(2));
        lambda  = pars(3);
        
        ns = [6 9 12 15];
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(4); p_aff9 = pars(5); p_aff12 = pars(6); p_aff15 = pars(7); 
            
            for i_n = 1:length(ns)
                filter_n = find(stim.N == ns(i_n));

                for i_filtn = 1:length(filter_n)
                    trial_idx = filter_n(i_filtn);
                    X_n{i_filtn} = stim.X{trial_idx};
                    Y_n{i_filtn} = stim.Y{trial_idx};
                end

                false_d(filter_n) = get_false_affiliation_bayesian_d(X_n, Y_n, eval(p_affs{i_n}));
            end
         else
             false_d = stim.bayesd;
        end
        %%%
        
        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), k, sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 2 % (A) 1 with 4 K's  --Bayes 
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            for i_n = 1:length(ns)
                filter_n = find(stim.N == ns(i_n));

                for i_filtn = 1:length(filter_n)
                    trial_idx = filter_n(i_filtn);
                    X_n{i_filtn} = stim.X{trial_idx};
                    Y_n{i_filtn} = stim.Y{trial_idx};
                end

                false_d(filter_n) = get_false_affiliation_bayesian_d(X_n, Y_n, eval(p_affs{i_n}));
            end
         else
             false_d = stim.bayesd;
        end
        %%%
        
        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 3 % (A) 1 with false sigma_s and sigma -- false Bayes
        k       = pars(1);
        sigma_d = exp(pars(2));
        lambda  = pars(3);
        sigma_s = exp(pars(4));
        sigma   = exp(pars(5));        
        
        false_d = get_false_bayesian_d(stim.X, stim.Y, sigma_s, sigma);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(false_d, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 4 % (B) Centroid (mean density) -- mu-conditioned, z-marginalized
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_meandist_mu(x_n, y_n, eval(p_affs{n_idx}));
            end
         else
             false_d = stim.MeanDist_Mu.LLR_mu;
        end
        %%%
        
        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 5 % (B) Max local density-- mu-conditioned, z-marginalized
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_localgaussgridmu(x_n, y_n, eval(p_affs{n_idx}));
            end
         else
             false_d = stim.localgaussgrid.LLR;
        end
        %%%
        
        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
       
        
    case 6 % (B) mu_PME -- mu-conditioned, z-marginalized
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_mu_PME(x_n, y_n, eval(p_affs{n_idx}));
            end
         else
             false_d = stim.PME_mu.LLR_mu;
        end
        %%%


        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 7 % Max marg likelihood (mostlikelyz)
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_mostlikelyz(x_n, y_n, eval(p_affs{n_idx}));
            end
         else
             false_d = stim.mostlikelyz.LLR;
        end
        %%%
        
        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d); %%%
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 8 % Agglomerative clustering
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);
        
        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_centroidadd(x_n, y_n, eval(p_affs{n_idx}));
            end
         else
             false_d = stim.centroidadd.LLR_z;
        end
        %%%

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d); %%%
        end
  
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 9 % (D) Brute max joint posterior
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
%             for trial_idx = 1:length(stim.X)
%                 n_idx = find(ns == stim.N(trial_idx));
%                 
%                 x_n = stim.X{trial_idx}';
%                 y_n = stim.Y{trial_idx}';
% 
%                 [~, false_d(trial_idx)] = func_brute_max(x_n, y_n, eval(p_affs{n_idx}));
%             end

             for trial_idx = 1:length(stim.X)
                 n_idx = find(ns == stim.N(trial_idx));
                 n_ = ns(n_idx);
                 precomputed_LLR_trial = precomputed_LLR.post{trial_idx};

                 false_d(trial_idx) = func_brute_max_precompute(n_, eval(p_affs{n_idx}), precomputed_LLR_trial); 
             end
            
         else
             false_d = stim.brutemax.posterior.LLR_z_mu;
        end
        %%%

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
 
     case 10 % (D) Brute max joint likelihood
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 

             for trial_idx = 1:length(stim.X)
                 n_idx = find(ns == stim.N(trial_idx));
                 n_ = ns(n_idx);
                 precomputed_LLR_trial = precomputed_LLR.like{trial_idx};
                 
                 false_d(trial_idx) = func_brute_max_precompute(n_, eval(p_affs{n_idx}), precomputed_LLR_trial); 
             end
             
         else
             false_d = stim.brutemax.likelihood.LLR_z_mu;
        end
        %%%

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
  
     case 11 % (D) Centroid add-one-- z-conditioned and mu-conditioned
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_centroidadd_z_mu(x_n, y_n, eval(p_affs{n_idx}));
            end
            
         else
             false_d = stim.centroidadd_z_mu.LLR_z_mu;
        end
        %%%

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 12 % (D) Centroid (Mean density) -- mu-conditioned, z-optimal threshold conditioned
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
        %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_maxz_meandistmu(x_n, y_n, eval(p_affs{n_idx}));
            end
            
         else
             false_d = stim.MeanDist_Mu.LLR_z_mu;
        end
        
        %%%
        
        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
      case 13 % (D) max local dens, mu-conditioned, z-opt conditioned
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};
        
         %%%
        false_d = nan(length(stim.X),1);
        
        if paffmodel == 1
            p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
            p_aff6 = pars(7); p_aff9 = pars(8); p_aff12 = pars(9); p_aff15 = pars(10); 
            
            
            for trial_idx = 1:length(stim.X)
                n_idx = find(ns == stim.N(trial_idx));
                
                x_n = stim.X{trial_idx}';
                y_n = stim.Y{trial_idx}';

                false_d(trial_idx) = func_maxz_localgaussgridmu(x_n, y_n, eval(p_affs{n_idx}));
            end
            
         else
             false_d = stim.localgaussgrid.LLR_z_mu;
        end
        %%%

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 14 % Fraction of points (z) below NNdist threshold (d)
        d6      = exp(pars(1));
        d9      = exp(pars(2));
        d12     = exp(pars(3));
        d15     = exp(pars(4));
        z6      = pars(5);
        z9      = pars(6);
        z12     = pars(7);
        z15     = pars(8);
        lambda  = pars(9);

        ns = [6 9 12 15];
        ds = {'d6' 'd9' 'd12' 'd15'};
        zs = {'z6' 'z9' 'z12' 'z15'};

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            d = eval(ds{i_n});
            z = eval(zs{i_n});
            
            for i_trial = 1:length(filter_n)
                trialidx = filter_n(i_trial);
                
                fractbelow = sum(stim.NNdist{trialidx}<d)./stim.N(trialidx);
                p_RF(trialidx) = lambda/2 + (1-lambda)*(fractbelow > z);
            end
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 15 % Local density -- Gaussian convolve
        k6      = exp(pars(1));
        k9      = exp(pars(2));
        k12     = exp(pars(3));
        k15     = exp(pars(4));
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(stim.localgaussgrid.max(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 16 % mean NN dist
        k6      = exp(pars(1));
        k9      = exp(pars(2));
        k12     = exp(pars(3));
        k15     = exp(pars(4));
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RNF(filter_n) = lambda/2 + (1-lambda) * normcdf(stim.meanNNdist(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RF = 1 - p_RNF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 17 % min dist
        k6      = exp(pars(1));
        k9      = exp(pars(2));
        k12     = exp(pars(3));
        k15     = exp(pars(4));
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RNF(filter_n) = lambda/2 + (1-lambda) * normcdf(stim.MinDist(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RF = 1 - p_RNF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 18 % local dens
        k6      = exp(pars(1));
        k9      = exp(pars(2));
        k12     = exp(pars(3));
        k15     = exp(pars(4));
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(stim.localgaussgrid.max(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 19 % global dens
        k6      = exp(pars(1));
        k9      = exp(pars(2));
        k12     = exp(pars(3));
        k15     = exp(pars(4));
        sigma_d = exp(pars(5));
        lambda  = pars(6);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RNF(filter_n) = lambda/2 + (1-lambda) * normcdf(stim.MeanDist(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RF = 1 - p_RNF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
end
end

                    