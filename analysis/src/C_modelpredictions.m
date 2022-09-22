function [NLL, p_RF, false_d] = C_modelpredictions(pars, modelidx, stim, data, paffmodel, precomputed_LLR)

% Fitting mode or model prediction mode?
nodata      = isempty(data);

switch modelidx
    case 1 % bayes 
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.bayesd, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 2 % global dens
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
    case 3 % local dens
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
     case 4 % min dist
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 5 % mean NN dist
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 6 % EM-based mu-conditioned simulated opt-k LLR
        % sim_optk.EM
        k6      = 2.4242;
        k9      = 2.7273;
        k12     = 3.0303;
        k15     = 3.0303;
        sigma_d = pars(1);
        lambda  = pars(2);

        ns = [6 9 12 15];
        ks = {'k6' 'k9' 'k12' 'k15'};

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(stim.EM.LLR(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
      case 7 % Model  7: EM-based mu-conditioned LLR where k is fit for all n's
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.EM.LLR, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
        
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
      case 8 % Model  8: Local dens-based LLR where k is fit for all n's
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.localgaussgrid.LLR, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
        
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
      case 9 % Model  9: Global mean dist then bayes (hist ratio)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.simLLR.MeanDist, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
      case 10 % Model  10: Min distance then bayes (hist ratio)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.simLLR.MinDist, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
      case 11 % Model  11: Mean NN distance then bayes (hist ratio)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.simLLR.meanNNdist, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
    case 12 % Model  12: (Max) Local density then bayes (hist ratio)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.simLLR.localgaussgrid, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 13 % Model  13: EM then bayes (hist ratio)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.simLLR.EM, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 14 % Model  14: (Sum) Local density -- Gaussian convolve
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
     
     case 15 % Model  15: Gaussian radius likelihood threshold z-conditioned, mean dist mu-conditioned
        k           = pars(1);
        sigma_d     = pars(2);
        lambda      = pars(3);
        like_thresh = pars(4);
        
        % Get z-vectors based on likelihood thresholding, then get LLR_z_mu
        for i_trial = 1:length(stim.X)
            x   = stim.X{i_trial};
            y   = stim.Y{i_trial};
            mu  = stim.MeanDist_Mu.trial(i_trial,:);
            z   = stim.MeanDist_Mu.C1_likelihoods{i_trial} > like_thresh;
            
            LLR_z_mu(i_trial) = get_LLR_z_mu(x,y,z,mu);
        end
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 16 % Model  16: Fraction of points (z) below NNdist threshold (d)
        d6      = pars(1);
        d9      = pars(2);
        d12     = pars(3);
        d15     = pars(4);
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
        
        
     case 17 % Model  17: Centroid-add-one (BHC)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.centroidadd.LLR_z, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end

        
     case 18 % Model  18: Centroid add-one-- z-conditioned and mu-conditioned
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.centroidadd_z_mu.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end

     case 19 % Model  19: Centroid add-one-- z-conditioned throughout, with FINAL z- and mu-conditioning
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.centroidadd.LLR_z_finalmu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end

     case 20 % Model  20: Centroid add-one-- z- and mu-conditioned throughout, with FINAL z-conditioning
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.centroidadd_z_mu.LLR_z_mu_finalmumarg, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 21 % Model  21: Hard EM LLR_z (z-conditioned, mu-marginalized)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.EM.hard.LLR_z, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 22 % Model  22: Hard EM LLR_z_mu (z-conditioned and mu-conditioned)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.EM.hard.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 23 % Model  23: Gaussian radius likelihood threshold z-conditioned, EM mu-conditioned
        k           = pars(1);
        sigma_d     = pars(2);
        lambda      = pars(3);
        like_thresh = pars(4);
        
        % Get z-vectors based on likelihood thresholding, then get LLR_z_mu
        for i_trial = 1:length(stim.X)
            x   = stim.X{i_trial};
            y   = stim.Y{i_trial};
            mu  = stim.EM.mu(i_trial,:);
            z   = stim.EM.C1_likelihoods{i_trial} > like_thresh;
            
            LLR_z_mu(i_trial) = get_LLR_z_mu(x,y,z,mu);
        end
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 24 % Model  24: Gaussian radius likelihood threshold z-conditioned based on EM mu, mu-marginalized
        k           = pars(1);
        sigma_d     = pars(2);
        lambda      = pars(3);
        like_thresh = pars(4);
        
        % Get z-vectors based on likelihood thresholding, then get LLR_z --
        % mu marg
        for i_trial = 1:length(stim.X)
            x   = stim.X{i_trial};
            y   = stim.Y{i_trial};
            z   = stim.EM.C1_likelihoods{i_trial} > like_thresh;
            
            LLR_z(i_trial) = get_LLR_z(x,y,z);
        end
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(LLR_z, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end   
        
     case 25 % Model  25: Gaussian radius likelihood threshold z-conditioned based on mean dist mu mu-marginalized (cf 15)
        k           = pars(1);
        sigma_d     = pars(2);
        lambda      = pars(3);
        like_thresh = pars(4);
        
        % Get z-vectors based on likelihood thresholding, then get LLR_z --
        % mu marg
        for i_trial = 1:length(stim.X)
            x   = stim.X{i_trial};
            y   = stim.Y{i_trial};
            z   = stim.MeanDist_Mu.C1_likelihoods{i_trial} > like_thresh;
            
            LLR_z(i_trial) = get_LLR_z(x,y,z);
        end
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(LLR_z, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end   
        
    case 26 % Model  26: Most likely z-- z estimated by maximizing marginal likelihood of z (marginal posterior is the same)
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.mostlikelyz.LLR, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end        

    case 27 % Model  27: Max marginal likelihood mu-conditioned, z-marginalized
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.maxmargmu.likelihood.LLR, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end        

     case 28 % Model  28: Max marginal posterior mu-conditioned, z-marginalized
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.maxmargmu.posterior.LLR, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end   
        
     case 29 % Model  29: Max local density-- mu-conditioned, z-marginalized
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.localgaussgrid.LLR, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 30 % Model  30: Mean density -- mu-conditioned, z-marginalized
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.MeanDist_Mu.LLR_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end    
        
     case 31 % Model  31: Soft EM -- mu-conditioned, z-marginalized
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.EM.LLR, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end 
        
     case 32 % Model  32: Max marginal likelihood mu-conditioned, z-optimal threshold conditioned 
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.maxmargmu.likelihood.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end      
        
     case 33 % Model  33: Max marginal posterior mu-conditioned, z-optimal threshold conditioned
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.maxmargmu.posterior.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end                
   
     case 34 % Model  34: Max local density -- mu-conditioned, z-optimal threshold conditioned
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.localgaussgrid.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 35 % Model  35: Mean density -- mu-conditioned, z-optimal threshold conditioned
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.MeanDist_Mu.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 36 % Model  36: Soft EM -- mu-conditioned, z-optimal threshold conditioned
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.EM.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 37 % Model  37: Brute max joint likelihood
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.brutemax.likelihood.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
        
     case 38 % Model  38: Brute max joint posterior
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(stim.brutemax.posterior.LLR_z_mu, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end       
        
      case 39 % Model  39: (C+) 17 with 4 K's
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 40 % Model  40: (C+) 26 with 4 K's
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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

     case 41 % Model  41: (D+) 32 with 4 K's % RETIRED
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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

                false_d(trial_idx) = func_maxz_maxmargmu_likelihood(x_n, y_n, eval(p_affs{n_idx}));
            end
            
         else
             false_d = stim.maxmargmu.likelihood.LLR_z_mu;
        end
        %%%
        
        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            p_RF(filter_n) = lambda/2 + (1-lambda) * normcdf(false_d(filter_n), eval(ks{i_n}), sigma_d);
        end
        p_RNF = 1 - p_RF;
        
        p_RF((p_RF==0)) = 1e-300;
        p_RNF((p_RNF==0)) = 1e-300;

        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end

      case 42 % Model  42: (D+) 33 with 4 K's % RETIRED
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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

                false_d(trial_idx) = func_maxz_maxmargmu_posterior(x_n, y_n, eval(p_affs{n_idx}));
            end
            
         else
             false_d = stim.maxmargmu.posterior.LLR_z_mu;
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

      case 43 % Model  43: (D+) 34 with 4 K's -- max local dens, mu-conditioned, z-opt conditioned
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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

      case 44 % Model  44: (D+) 35 with 4 K's Mean density -- mu-conditioned, z-optimal threshold conditioned
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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

      case 45 % Model  45: (D+) 18 with 4 K's --Centroid add-one-- z-conditioned and mu-conditioned
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 46 % Model  46: (D+) 37 with 4 K's --Brute max joint likelihood
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 47 % Model  47: (D+) 38 with 4 K's --Brute max joint posterior
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 48 % Model  48: (B+) 27 with 4 K's --Max marginal likelihood-- mu-conditioned, z-marginalized
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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

                false_d(trial_idx) = func_maxmargmu_like(x_n, y_n, eval(p_affs{n_idx}));
            end
         else
             false_d = stim.maxmargmu.likelihood.LLR;
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
        
      case 49 % Model  49: (B+) new model -- mu_PME -- mu-conditioned, z-marginalized
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 50 % Model  50: (B+) 29 with 4 K's --Max local density-- mu-conditioned, z-marginalized
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 51 % Model  51: (B+) 30 with 4 K's --Mean density -- mu-conditioned, z-marginalized
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 52 % Model  52: (A+) 1 with 4 K's  --Bayes 
        k6      = pars(1);
        k9      = pars(2);
        k12     = pars(3);
        k15     = pars(4);
        sigma_d = pars(5);
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
        
      case 53 % Model  53: (A) 1 with false sigma_s and sigma -- false Bayes
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        sigma_s = pars(4);
        sigma   = pars(5);        
        
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
        
      case 54 % Model  54: (A) 1 with false affiliation belief -- false Bayes
        k       = pars(1);
        sigma_d = pars(2);
        lambda  = pars(3);
        p_aff6  = pars(4);
        p_aff9  = pars(5);
        p_aff12 = pars(6);
        p_aff15 = pars(7);

        ns = [6 9 12 15];
        p_affs = {'p_aff6' 'p_aff9' 'p_aff12' 'p_aff15'};
        false_d = nan(length(stim.X),1);

        for i_n = 1:length(ns)
            filter_n = find(stim.N == ns(i_n));
            
            for i_filtn = 1:length(filter_n)
                trial_idx = filter_n(i_filtn);
                X_n{i_filtn} = stim.X{trial_idx};
                Y_n{i_filtn} = stim.Y{trial_idx};
            end
            
            false_d(filter_n) = get_false_affiliation_bayesian_d(X_n, Y_n, eval(p_affs{i_n}));
        end
        
        p_RF    = lambda/2 + (1-lambda) * normcdf(false_d, k, sigma_d); % cdf evaluated at d, with mean k and sigma_d
        p_RNF   = 1 - p_RF;

        p_RF(find(p_RF==0)) = 1e-300;
        p_RNF(find(p_RNF==0)) = 1e-300;
         
        if nodata
            NLL = [];
        else
            NLL = -(sum(log(p_RF(find(data.Resp_feeder)))) + sum(log(p_RNF(find(~data.Resp_feeder)))));
        end
end

end
