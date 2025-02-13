% #########################################################################
%
% Simultaneous state estimation and control for nonlinear systems under
% bounded disturbances.
%
% #########################################################################
%
% #########################################################################


function sim = vanDerPol_mulY()
    clear all; clc;% close all;
    import casadi.*
    sim = [];
    %
    sim = init(sim);    
    sim = gen_constraints(sim);
    %
    for k=sim.config.NeNe
        sim = set_Ne(sim, k);
        for m=4%:2:6%&10
            sim = gen_noise_variance(sim, m, 2, 0.25, 0);
            for l=1:length(sim.config.densities)
                sim = set_density_outliers(sim, sim.config.densities(l));
                for j=1:sim.config.numBatchs                    
                    sim = gen_c2(sim);                      % Used with mheLoss
                    for i=1:sim.config.numTrials
                        sim = gen_x0(sim);
                        sim = gen_noise_sequences(sim);     % place outside in ordert to generate always the same state trajectory
                        sim = gen_state_trajectory(sim);
                        %
                        sim = gen_x0bar(sim);
                        sim = build_estimators(sim);
                        sim = set_init_condition_estimators(sim);
                        sim = estimate(sim);
                        sim = append_to_plot(sim);
                        %
                        a = [k, m, l, j, i]
                    end
                    sim = compute_errors(sim);
                end
            end
        end
    end
    %
    sim = plt(sim);
end

function sim = build_setup(sim)
    % SIMULATION PARAMETERS ===================================================
    sim.config.t0                   = 0;                                                 % initial time of sym. [Seg]
    sim.config.tf                   = 10;                                                % final time of sym. [Seg]
    sim.config.Ts                   = 0.025;                                             % sampling period, 0=>discrete-time system
    sim.config.N                    = (sim.config.tf-sim.config.t0)/sim.config.Ts-1;     % number of steps of sym.
    sim.config.same_seed            = false;
    sim.config.Ne                   = 10;
    sim.config.NeNe                 = 10;%[5 10 15 20];
    % Noise amplitudes __________________________________________________________
    sim.data.processDistVariance    = [0.01; 0.015];                                     % process noise amplitude
    sim.data.measNoiseVariance      = [0.25; 0.25; 0.45; 0.45];%; 0.35; 0.35; 0.35; 0.35; 0.35; 0.35];             % measurement noise amplitude: nv x numSensors     
    sim.data.measNoiseMean          = [0; 0; 0; 0];%; 0; 0; 0; 1; 0; 0];             % measurement noise amplitude: nv x numSensors     
    sim.data.outlier                = true;
    sim.config.maxChannelsWithOutliers = 2;                                          % [1,length(sim.data.measNoiseVariance)]; NMHE takes only the lats channel measurements, include outliers in this chanel for properly compare all methods
    sim.config.amplOutlier          = 10;   % x times the value of sim.data.measNoiseVariance(sim.config.maxChannelsWithOutliers)
    sim.config.densityOutliers      = 0.35;  % value in [0, 1]: 0: no outliers at channel "sim.config.maxChannelsWithOutliers", 1: every sample in channel "sim.config.maxChannelsWithOutliers" with oputliers
    sim.config.densities            = 0.3;%:0.1:0.9;
    sim.config.noiseDist            = 'gauss';% 'uniform;' PŔOBAR CON UNA DISTRIBUCION BIMODAL
    %
    sim.config.Max_err_x0           = 10;                                                % deviation of the initial guess
    %
    sim.time.total_sim              = 0;
    %
    sim.config.numTrials            = 1;
    sim.config.numBatchs            = 1;
    if sim.config.numTrials > 1
        sim.config.same_seed = false;
    end
    %
    sim.config.numSensors           = 5;
    %
    sim.config.fractionTs           = 10;                        % percentual value: [0, 99]
    sim.config.delayY               = true;
    %    
end

function sim = set_fraction_of_Ts(sim, fraction)
    sim.config.fractionTs = fraction;
end

function sim = set_density_outliers(sim, density)
    sim.config.densityOutliers = density;
end

function sim = set_Ne(sim, Ne)
    sim.config.Ne = Ne;
end

function sim = gen_noise_variance(sim, numChannel, numChannelOutliers, dev, mean)       % numChannel = nx*num_sensors
    sim.data.measNoiseVariance          = repmat(dev, numChannel, 1);
    sim.data.measNoiseMean              = repmat(mean, numChannel, 1);
    sim.config.maxChannelsWithOutliers  = numChannelOutliers;
    %
    sim.config.numSensors               = numChannel/sim.system.nv;
    %
    sim.data.ysim1          = zeros(sim.system.ny*sim.config.numSensors,sim.config.N);
end

function sim = gen_c2(sim)
%     c2                  = rand*0.1;
%     c2                  = rand*0.01;
    c2                  = rand*0.001;
    % c2                  = rand*0.0001;
    sim.data.mheLoss.c2 = [sim.data.mheLoss.c2; c2];
end

function sim = gen_dynamic(sim)
    % System Parameters _______________________________________________________
    sim.system.nx      = 2;                            % vector state dimension
    sim.system.nw      = 2;                            % process noise input dimension
    sim.system.nd      = 2;                            % model mismatch disturbance
    sim.system.ny      = 2;                            % output dimension
    sim.system.nu      = 1;                            % input dimension
    sim.system.nv      = sim.system.ny;
    % Casadi variabes ----------------------------------------------------
    sim.dynamic.x      = casadi.MX.sym('x',sim.system.nx);
    sim.dynamic.w      = casadi.MX.sym('w',sim.system.nw);
    sim.dynamic.u      = casadi.MX.sym('u',sim.system.nu);
    % System's parameters ------------------------------------------------
    sim.system.epsilon = 0;
    sim.system.freq    = 3;
    % Dynamic of the system ----------------------------------------------
    sim.dynamic.f_rhs  = [sim.system.epsilon*(1 - sim.dynamic.x(2)^2)*sim.dynamic.x(1) - sim.system.freq*sim.dynamic.x(2) + sim.dynamic.u; sim.system.freq*sim.dynamic.x(1)];
    sim.dynamic.f      = casadi.Function('f', {sim.dynamic.x,sim.dynamic.u}, {sim.dynamic.f_rhs});  
    sim.dynamic.fToLin = sim.dynamic.f;    
    % Output of the system ------------------------------------------------
    sim.dynamic.h_rhs  = [sim.dynamic.x(1); sim.dynamic.x(2)];
    sim.dynamic.h      = casadi.Function('h', {sim.dynamic.x}, {sim.dynamic.h_rhs});  
    % compute jacobians ---------------------------------------------------
    sim.dynamic.jac_fx = casadi.Function('Func_A', {sim.dynamic.x, sim.dynamic.u}, {sim.dynamic.f_rhs.jacobian(sim.dynamic.x)}, {'x', 'u'}, {'A'});
    sim.dynamic.jac_fu = casadi.Function('Func_B', {sim.dynamic.x, sim.dynamic.u}, {sim.dynamic.f_rhs.jacobian(sim.dynamic.u)}, {'x', 'u'}, {'B'});
    sim.dynamic.jac_hx = casadi.Function('Func_C', {sim.dynamic.x}, {sim.dynamic.h_rhs.jacobian(sim.dynamic.x)}, {'x'}, {'C'});
    sim.dynamic.jac_fhx = casadi.Function('Func_C', {sim.dynamic.x}, {sim.dynamic.f_rhs.jacobian(sim.dynamic.x)}, {'x'}, {'C'});
    % RK4 -----------------------------------------------------------------
    k1          = sim.dynamic.f(sim.dynamic.x,sim.dynamic.u);
    k2          = sim.dynamic.f(sim.dynamic.x + sim.config.Ts / 2 * k1, sim.dynamic.u);
    k3          = sim.dynamic.f(sim.dynamic.x + sim.config.Ts / 2 * k2, sim.dynamic.u);
    k4          = sim.dynamic.f(sim.dynamic.x + sim.config.Ts * k3, sim.dynamic.u);
    x_rk4       = sim.dynamic.x + sim.config.Ts / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    sim.dynamic.F = casadi.Function('F', {sim.dynamic.x, sim.dynamic.u}, {x_rk4});
    %
    k1          = sim.dynamic.f(sim.dynamic.x,sim.dynamic.u);
    k2          = sim.dynamic.f(sim.dynamic.x + ((sim.config.Ts * (100-sim.config.fractionTs)/100)) / 2 * k1, sim.dynamic.u);
    k3          = sim.dynamic.f(sim.dynamic.x + ((sim.config.Ts * (100-sim.config.fractionTs)/100)) / 2 * k2, sim.dynamic.u);
    k4          = sim.dynamic.f(sim.dynamic.x + ((sim.config.Ts * (100-sim.config.fractionTs)/100)) * k3, sim.dynamic.u);
    x_rk4_MulYDelay   = sim.dynamic.x + ((sim.config.Ts * (100-sim.config.fractionTs)/100)) / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    sim.dynamic.FMulYDelay = casadi.Function('FMulYDelay', {sim.dynamic.x, sim.dynamic.u}, {x_rk4_MulYDelay});
    %
    k1          = sim.dynamic.f(sim.dynamic.x,sim.dynamic.u);
    k2          = sim.dynamic.f(sim.dynamic.x + ((sim.config.Ts * sim.config.fractionTs/100)/(sim.config.numSensors-1)) / 2 * k1, sim.dynamic.u);
    k3          = sim.dynamic.f(sim.dynamic.x + ((sim.config.Ts * sim.config.fractionTs/100)/(sim.config.numSensors-1)) / 2 * k2, sim.dynamic.u);
    k4          = sim.dynamic.f(sim.dynamic.x + ((sim.config.Ts * sim.config.fractionTs/100)/(sim.config.numSensors-1)) * k3, sim.dynamic.u);
    x_rk4_MulY  = sim.dynamic.x + ((sim.config.Ts * sim.config.fractionTs/100)/(sim.config.numSensors-1)) / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    sim.dynamic.FMulY = casadi.Function('FMulY', {sim.dynamic.x, sim.dynamic.u}, {x_rk4_MulY});
    % Used only with yhe bounds calculator --------------------------------
    sim.dynamic_anon.f_anon  = @(t,x,w,u)[sim.system.epsilon*(1 - x(2)^2)*x(1) - sim.system.freq*x(2) + w(1) + u; sim.system.freq*x(1) + w(2)];
    sim.dynamic.anon.h_anon  = @(x,v)[x(1) + x(2); x(1) - x(2)];
end

function sim = init(sim)
    sim = build_setup(sim);
    sim = gen_dynamic(sim);    
    sim = reserve_memory(sim);
end

function sim = build_estimators(sim)
    sim = build_nlmhe(sim);
    sim = build_nlmheAvg(sim);
    sim = build_nlmheMulY(sim);
    sim = build_nlmheConvexY(sim);
    sim = build_ekf_joint(sim);
    sim = build_nlmheLoss(sim);
    sim = build_nlmheAlessandri(sim);
end

function sim = append_to_plot(sim)
    sim.data.ysim_plt           = [sim.data.ysim_plt; sim.data.ysim1];
    %
    sim.data.xsim_plt           = [sim.data.xsim_plt;sim.data.xsim1];
    sim.data.xestmhe_plt        = [sim.data.xestmhe_plt;sim.data.xestmhe];
    sim.data.xestmheAvg_plt     = [sim.data.xestmheAvg_plt;sim.data.xestmheAvg];
    sim.data.xestmheMulY_plt    = [sim.data.xestmheMulY_plt;sim.data.xestmheMulY];
    sim.data.xestmheConvexY_plt = [sim.data.xestmheConvexY_plt;sim.data.xestmheConvexY];
    sim.data.xestekf_plt        = [sim.data.xestekf_plt;sim.data.xestjoint];
    sim.data.xestmheLoss_plt    = [sim.data.xestmheLoss_plt;sim.data.xestmheLoss];
    sim.data.xestmheAle_plt     = [sim.data.xestmheAle_plt;sim.data.xestmheAle];
    %
    sim.data.timemhe_plt        = [sim.data.timemhe_plt;sim.data.timemhe];
    sim.data.timemheAvg_plt     = [sim.data.timemheAvg_plt;sim.data.timemheAvg];
    sim.data.timemheMulY_plt    = [sim.data.timemheMulY_plt;sim.data.timemheMulY];
    sim.data.timemheConvexY_plt = [sim.data.timemheConvexY_plt;sim.data.timemheConvexY];
    sim.data.timeekf_plt        = [sim.data.timeekf_plt;sim.data.timeJoint];
    sim.data.timemheLoss_plt    = [sim.data.timemheLoss_plt;sim.data.timemheLoss];
    sim.data.timemheAle_plt     = [sim.data.timemheAle_plt;sim.data.timemheAle];
end

function sim = set_init_condition_estimators(sim)
%     set_x0bar(sim.ekf.ekf_joint, sim.estimator.x0bar);
%     set_alpha0bar(sim.ekf.ekf_joint, sim.polytope.alpha0bar);
    set_x0bar(sim.mhe.nlmhe, sim.estimator.x0bar);
    %
    set_x0bar(sim.mheAvg.nlmheAvg, sim.estimator.x0bar);
    %
    set_x0bar(sim.mheMulY.nlmheMulY, sim.estimator.x0bar);
    %
    set_x0bar(sim.mheConvexY.nlmheConvexY, sim.estimator.x0bar);
    %
    set_x0bar(sim.mheLoss.nlmheLoss, sim.estimator.x0bar);
    %
    set_x0bar(sim.mhe.nlmheAle, sim.estimator.x0bar);
end

function sim = reserve_memory(sim)
    % Number of simulation to carry-on ____________________________________
    sim.data.Jmhe           = casadi.DM.zeros(1,sim.config.N);
    sim.data.JmheAvg        = casadi.DM.zeros(1,sim.config.N);
    sim.data.JmheMulY       = casadi.DM.zeros(1,sim.config.N);
    sim.data.JmheConvexY    = casadi.DM.zeros(1,sim.config.N);
    sim.data.JMMmhe         = [];
    sim.data.JmheLoss       = casadi.DM.zeros(1,sim.config.N);
    sim.data.arrCost        = [];
    % Variables: states, outputs, etc. ____________________________________
    sim.data.xsim1          = zeros(sim.system.nx,sim.config.N+1);
    sim.data.ysim1          = zeros(sim.system.ny*sim.config.numSensors,sim.config.N);
    sim.data.Uk1            = zeros(1,sim.config.N);%0.1*idinput([20,1,(sim.config.N+1)/5])';
    sim.data.xsim1Extra     = zeros(sim.system.nx,(sim.config.N+1)*sim.config.numSensors);
    %
    sim.data.estErrs            = [];
    sim.data.estErrsGlb         = [];
    %
    sim.data.ysim_plt           = [];
    sim.data.xsim_plt           = [];
    sim.data.xestmhe_plt        = [];
    sim.data.xestmheAvg_plt     = [];
    sim.data.xestmheMulY_plt    = [];
    sim.data.xestmheConvexY_plt = [];
    sim.data.xestekf_plt        = [];
    sim.data.xestmheLoss_plt    = [];
    sim.data.xestmheAle_plt     = [];
    %
    sim.data.mheLoss.alpha_w = [];
    sim.data.mheLoss.alpha_v = [];
    %
    sim.data.xestmhe         = NaN(sim.system.nx,sim.config.N);
    sim.data.xestmheAvg      = NaN(sim.system.nx,sim.config.N);
    sim.data.xestmheMulY     = NaN(sim.system.nx,sim.config.N);
    sim.data.xestmheConvexY  = NaN(sim.system.nx,sim.config.N);
    sim.data.xestjoint       = NaN(sim.system.nx,sim.config.N);
    sim.data.xestmheLoss     = NaN(sim.system.nx,sim.config.N);
    sim.data.xestmheAle      = NaN(sim.system.nx,sim.config.N);
    %
    sim.data.timemhe_plt        = [];
    sim.data.timemheAvg_plt     = [];
    sim.data.timemheMulY_plt    = [];
    sim.data.timemheConvexY_plt = [];
    sim.data.timeekf_plt        = [];
    sim.data.timemheLoss_plt    = [];
    sim.data.timemheAle_plt     = [];
    %
    sim.data.timemhe         = NaN(1,sim.config.N);
    sim.data.timemheAvg      = NaN(1,sim.config.N);
    sim.data.timemheMulY     = NaN(1,sim.config.N);
    sim.data.timemheConvexY  = NaN(1,sim.config.N);
    sim.data.timeJoint       = NaN(1,sim.config.N);
    sim.data.timemheLoss     = NaN(1,sim.config.N);
    sim.data.timemheAle      = NaN(1,sim.config.N);
    %
    sim.data.Pinv_mhe        = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    sim.data.Pinv_mheAvg     = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    sim.data.Pinv_mheMulY    = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    sim.data.Pinv_mheConvexY = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    sim.data.Pinv_mheLoss    = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    %
    sim.data.max_eig_P_mhe   = NaN(1,sim.config.N);
    sim.data.max_eig_P_mheAvg = NaN(1,sim.config.N);
    sim.data.max_eig_P_mheMulY = NaN(1,sim.config.N);
    sim.data.max_eig_P_mheConvexY = NaN(1,sim.config.N);
    sim.data.max_eig_P_mheLoss = NaN(1,sim.config.N);
    %
    sim.data.min_eig_P_mhe   = NaN(1,sim.config.N);
    sim.data.min_eig_P_mheAvg = NaN(1,sim.config.N);
    sim.data.min_eig_P_mheMulY = NaN(1,sim.config.N);
    sim.data.min_eig_P_mheLoss = NaN(1,sim.config.N);sim.data.min_eig_P_mheConvexY = NaN(1,sim.config.N);
    %
    sim.data.mheLoss.c2         = [];
    sim.data.mheLoss.c          = [];
end

function sim = gen_constraints(sim)
    % Constraints _________________________________________________
    sim.estimator.constraints.Xe      = [];%[-5 5;-5 5];
    sim.estimator.constraints.We      = [];%3.*[-sim.data.processDistVariance sim.data.processDistVariance];
    sim.estimator.constraints.Ve      = [];%3.*[-sim.data.measNoiseVariance sim.data.measNoiseVariance];
    sim.estimator.constraints.Uc      = [];
    % only for the muliple model mhe    
end

function sim = gen_x0(sim)
    sim.system.x0       =  [0.58; -1.186];%&diag([3 -3]) * randn(2,1);% [3;0] + 2.*randn(2,1); [0.58; -1.186]
    sim.data.xsim1(:,1) = sim.system.x0;
end

function sim = gen_x0bar(sim)
    sim.estimator.x0bar = sim.system.x0 + (sim.config.Max_err_x0.*randn(sim.system.nx,1));
end

function sim = build_nlmhe(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mhe.l_x     = 1e6;
    sim.mhe.P       = sim.mhe.l_x.*casadi.DM.eye(sim.system.nx);        % init. val. of the arrival-cost weight matrix
    sim.mhe.Qx      = diag([15 15]);                                    % states stage-cost
    sim.mhe.Rx      = diag([10 10]);                                    % measurement noise stage-cost
    sim.mhe.RuP     = 1e6.*eye(sim.system.nu);                          % measurement noise stage-cost
    sim.mhe.Tsmhe   = sim.config.Ts;                                    % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mhe.nlmhe   = nlmheCasadi(sim.config.Ne,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mhe.Qx,sim.mhe.Rx,sim.mhe.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mhe.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
    setSigma(sim.mhe.nlmhe,10*max(max(sim.data.measNoiseVariance),1e-3));
    setC(sim.mhe.nlmhe,2*sim.mhe.l_x);
end

function sim = build_nlmheAlessandri(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheAle.l_x     = 1e6;
    sim.mheAle.P       = sim.mhe.l_x.*casadi.DM.eye(sim.system.nx);        % init. val. of the arrival-cost weight matrix
    sim.mheAle.Qx      = diag([15 15]);                                    % states stage-cost
    sim.mheAle.Rx      = diag([10 10]);                                    % measurement noise stage-cost
    sim.mheAle.Tsmhe   = sim.config.Ts;                                    % If 0, mhe does not integrate the systems
    sim.mheAle.beta    = 1;
    % *************************************************************
    sim.mhe.nlmheAle   = nlmheAlessandri(sim.config.Ne,0,sim.mheAle.beta,sim.dynamic.F,sim.dynamic.h,...
        sim.mheAle.Qx,sim.mheAle.Rx,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheAle.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts);

    set_x0bar(sim.mhe.nlmheAle,sim.estimator.x0bar);
    setJacobians(sim.mhe.nlmheAle,sim.dynamic.jac_fx,sim.dynamic.jac_fu,sim.dynamic.jac_hx);
%     obj = nlmheAlessandri(N,beta,x,w,u,f_rhs,h_rhs,Q,R,nx,nu,ny,nw,P,XLUBounds,WLUBounds,VLUBounds,ULUBounds,x0bar,Ts)
end

function sim = build_nlmheAvg(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheAvg.l_x     = 1e6;
    sim.mheAvg.P       = sim.mheAvg.l_x.*casadi.DM.eye(sim.system.nx);  % init. val. of the arrival-cost weight matrix
    sim.mheAvg.Qx      = diag([15 15]);                                 % states stage-cost
    sim.mheAvg.Rx      = diag([10 10]);                                 % measurement noise stage-cost
    sim.mheAvg.RuP     = 1e6.*eye(sim.system.nu);                       % measurement noise stage-cost
    sim.mheAvg.Tsmhe   = sim.config.Ts;                                 % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mheAvg.nlmheAvg   = nlmheCasadi(sim.config.Ne,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheAvg.Qx,sim.mheAvg.Rx,sim.mheAvg.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheAvg.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
    setSigma(sim.mheAvg.nlmheAvg,10*max(max(sim.data.measNoiseVariance),1e-3));
    setC(sim.mheAvg.nlmheAvg,2*sim.mheAvg.l_x);
end

function sim = build_nlmheMulY(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheMulY.l_x     = 1e6;
    sim.mheMulY.P       = sim.mheMulY.l_x.*casadi.DM.eye(sim.system.nx);    % init. val. of the arrival-cost weight matrix
    sim.mheMulY.Qx      = diag([15 15]);                                    % states stage-cost
    sim.mheMulY.Rx      = diag([10 10]);                                    % measurement noise stage-cost
    sim.mheMulY.RuP     = 1e6.*eye(sim.system.nu);                          % measurement noise stage-cost
    sim.mheMulY.Tsmhe   = sim.config.Ts;                                    % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mheMulY.nlmheMulY = nlmheMulY(sim.config.Ne,sim.config.numSensors,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheMulY.Qx,sim.mheMulY.Rx,sim.mheMulY.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheMulY.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
    setSigma(sim.mheMulY.nlmheMulY,10*max(max(sim.data.measNoiseVariance),1e-3));
    setC(sim.mheMulY.nlmheMulY,2*sim.mheMulY.l_x);
end

function sim = build_nlmheConvexY(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheConvexY.l_x     = 1e10;
    sim.mheConvexY.P       = sim.mheConvexY.l_x.*casadi.DM.eye(sim.system.nx); % init. val. of the arrival-cost weight matrix
    sim.mheConvexY.Qx      = diag([15 15]);                                 % states stage-cost
    sim.mheConvexY.Rx      = diag([10 10]);                                 % measurement noise stage-cost
    sim.mheConvexY.RuP     = 1e6.*eye(sim.system.nu);                       % measurement noise stage-cost
    sim.mheConvexY.Tsmhe   = sim.config.Ts;                                 % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mheConvexY.nlmheConvexY = nlmheConvexY(sim.config.Ne,sim.config.numSensors,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheConvexY.Qx,sim.mheConvexY.Rx,sim.mheConvexY.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheConvexY.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
    setSigma(sim.mheConvexY.nlmheConvexY,10*max(max(sim.data.measNoiseVariance),1e-3));
    setC(sim.mheConvexY.nlmheConvexY,2*sim.mheConvexY.l_x);
end

function sim = build_ekf_joint(sim)
%     sim.ekf.Qekf        = blkdiag(zeros(sim.system.nx),100.*eye(sim.mhe.mm.Nm));
%     sim.ekf.Rekf        = 0.01;
%     sim.ekf.Phat        = 1e-6.*blkdiag(eye(sim.system.nx),eye(sim.mhe.mm.Nm));
%     sim.ekf.ekf_joint   = kalman_x_alpha(sim.mhe.mm.Nm, sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.system.nv,...
%         sim.ekf.Phat,sim.mhe.mm.MMAdt,sim.mhe.mm.MMBdt,sim.mhe.mm.MMCdt,sim.mhe.mm.MMaidt,sim.estimator.x0bar,sim.polytope.alpha0,sim.ekf.Qekf,sim.ekf.Rekf);
end

function sim = build_nlmheLoss(sim)
    % ======================= ESTIMATOR ===========================         
%     sim.mheLoss.l_x     = 1e6;
%     sim.mheLoss.P       = sim.mhe.l_x.*casadi.DM.eye(sim.system.nx);        % init. val. of the arrival-cost weight matrix
%     sim.mheLoss.Qx      = diag([1 1]);                                    % states stage-cost
%     sim.mheLoss.Rx      = diag([1 1]);                                    % measurement noise stage-cost
%     sim.mheLoss.RuP     = 1e6.*eye(sim.system.nu);                          % measurement noise stage-cost
%     sim.mheLoss.Tsmhe   = sim.config.Ts;                                    % If 0, mhe does not integrate the systems
%     sim.mheLoss.c       = sim.data.mheLoss.c2(end);
%     % *************************************************************
% % nlmheCasadi_lossFunctions(N,x,w,u,f_rhs,h_rhs,cw,Kalphaw,cv,Kalphav,Ru,nx,nu,ny,nw,P,XLUBounds,WLUBounds,VLUBounds,ULUBounds,x0bar,integ,arrival_cost,arg_opt)    
%     sim.mheLoss.nlmheLoss   = nlmheCasadi_lossFunctions(sim.config.Ne,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
%         sim.mheLoss.c,sim.mheLoss.Qx,sim.mheLoss.Rx,sim.mhe.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mhe.P,sim.estimator.constraints.Xe,...
%         sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
%     setSigma(sim.mheLoss.nlmheLoss,10*max(max(sim.data.measNoiseVariance),1e-3));
%     setC(sim.mheLoss.nlmheLoss,2*sim.mheLoss.l_x);
    %
    dims        = {};
    dims.nq     = sim.system.nx;
    dims.nu     = sim.system.nu;
    dims.ny     = sim.system.ny;
    boxConst    = [];
    sim.mheLoss.nlmheLoss  = mheOptiLoss(sim.config.Ne, 0, sim.dynamic.F, sim.dynamic.h, dims, boxConst, sim.config.Ts);

    setSigmaP(sim.mheLoss.nlmheLoss,1.50); % 0.03
    setCP(sim.mheLoss.nlmheLoss,1e9);

    setJacobians(sim.mheLoss.nlmheLoss,sim.dynamic.jac_fx,sim.dynamic.jac_fu,sim.dynamic.jac_hx);
end

function sim = gen_noise_sequences(sim)    
    % Generate sequences of noises ________________________________
    if sim.config.same_seed
        sim.config_seed.s       = rng;              % seed number
        sim.config.seed.ran     = 134967344;        % set seed for rand generator
        rng(sim.config.seed.ran, 'twister');
    end
    sim.data.wc      = diag(sim.data.processDistVariance) * randn(sim.system.nw,sim.config.N);
    if strcmp(sim.config.noiseDist,'gauss')
        sim.data.vc      = diag(sim.data.measNoiseVariance) * randn(sim.system.nv*sim.config.numSensors,sim.config.N) + sim.data.measNoiseMean;
    elseif strcmp(sim.config.noiseDist,'uniform')
        sim.data.vc      = diag(sim.data.measNoiseVariance) * (2.*(rand(sim.system.nv*sim.config.numSensors,sim.config.N)-0.5)) + sim.data.measNoiseMean;
    end
    if sim.data.outlier
        posOutliers     = randperm(sim.config.N,ceil(sim.config.densityOutliers*sim.config.N));
        for i=1:length(posOutliers)
            channsOutliers  = randperm(length(sim.data.measNoiseVariance),sim.config.maxChannelsWithOutliers);
            sim.data.vc(channsOutliers,posOutliers(i)) = sim.config.amplOutlier * sim.data.vc(channsOutliers,posOutliers(i));
        end
    end
end

function sim = gen_state_trajectory(sim)    
    for i=1:sim.config.N
        % States __________________________________________________
        sim.data.xsim1(:,i+1) = full(sim.dynamic.F(sim.data.xsim1(:,i),sim.data.Uk1(1,i)))+sim.data.wc(:,i); 
        % Measurements ____________________________________________
        if sim.config.delayY
            if i==1
                for j=1:sim.config.numSensors
                    sim.data.ysim1((j-1)*sim.system.ny+1:j*sim.system.ny,i) = full(sim.dynamic.h(sim.data.xsim1(:,i))) + sim.data.vc((j-1)*sim.system.ny+1:j*sim.system.ny,i);
                    sim.data.xsim1Extra(:,(i-1)*sim.config.numSensors + j) = full(sim.data.xsim1(:,i));
                end
            else
                for j=1:sim.config.numSensors
                    if j==1
                        x = sim.dynamic.FMulYDelay(sim.data.xsim1(:,i-1),sim.data.Uk1(1,i-1))+sim.data.wc(:,i-1);
                    else
                        x = full(sim.dynamic.FMulY(x,sim.data.Uk1(1,i)))+sim.data.wc(:,i);
                    end
                    sim.data.ysim1((j-1)*sim.system.ny+1:j*sim.system.ny,i) = full(sim.dynamic.h(x)) + sim.data.vc((j-1)*sim.system.ny+1:j*sim.system.ny,i);
                    sim.data.xsim1Extra(:,(i-1)*sim.config.numSensors + j) = full(x);
                end
            end
        else
            for j=1:sim.config.numSensors
                sim.data.ysim1((j-1)*sim.system.ny+1:j*sim.system.ny,i) = full(sim.dynamic.h(sim.data.xsim1(:,i))) + sim.data.vc((j-1)*sim.system.ny+1:j*sim.system.ny,i);
            end
        end
    end
    %
    sim.data.xsim1 = sim.data.xsim1(:,1:sim.config.N);
    % plt
    xsim1                   = zeros(2,length(sim.data.xsim1)*sim.config.numSensors*floor(100/sim.config.fractionTs));    
    xsim1(:,sim.config.numSensors*floor(100/sim.config.fractionTs):sim.config.numSensors*floor(100/sim.config.fractionTs):end)    = sim.data.xsim1;    

    xsim1Extra              = zeros(2,length(sim.data.xsim1Extra)*floor(100/sim.config.fractionTs));    
    for i=sim.config.numSensors*floor(100/sim.config.fractionTs):sim.config.numSensors*floor(100/sim.config.fractionTs):length(xsim1Extra)
        xsim1Extra(:,i-(sim.config.numSensors-1):i) = sim.data.xsim1Extra(:,i/floor(100/sim.config.fractionTs)-(sim.config.numSensors-1):i/floor(100/sim.config.fractionTs));
    end

    ysim1Extra                   = zeros(2,length(sim.data.ysim1)*sim.config.numSensors*floor(100/sim.config.fractionTs));
    ysim1ExtraAux = reshape(sim.data.ysim1,2,sim.config.numSensors*length(sim.data.ysim1));
    for i=sim.config.numSensors*floor(100/sim.config.fractionTs):sim.config.numSensors*floor(100/sim.config.fractionTs):length(ysim1Extra)
        ysim1Extra(:,i-(sim.config.numSensors-1):i) = ysim1ExtraAux(:,i/floor(100/sim.config.fractionTs)-(sim.config.numSensors-1):i/floor(100/sim.config.fractionTs));
    end
    
    indx0                   = find(xsim1==0);
    xsim1(indx0)            = NaN;
    indx0                   = find(xsim1Extra==0);
    xsim1Extra(indx0)       = NaN;
    indx0                   = find(ysim1Extra==0);
    ysim1Extra(indx0)       = NaN;

    % figure;
    % subplot(2,1,1); hold on; grid on;
    % stem(ysim1Extra(1,:),'co','MarkerSize',4,'MarkerFaceColor','c');
    % stem(xsim1(1,:),'mo','MarkerSize',7,'MarkerFaceColor','m','LineWidth',1.25);
    % stem(xsim1Extra(1,:),'bo','MarkerSize',5,'MarkerFaceColor','b');    
    % 
    % subplot(2,1,2); hold on; grid on;
    % stem(ysim1Extra(2,:),'co','MarkerSize',4,'MarkerFaceColor','c');
    % stem(xsim1(2,:),'mo','MarkerSize',7,'MarkerFaceColor','m','LineWidth',1.25);
    % stem(xsim1Extra(2,:),'bo','MarkerSize',5,'MarkerFaceColor','b');    

end

function sim = estimate(sim)  
    for imhe=1:sim.config.N
        tic;
        if imhe < sim.config.Ne+1
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(end-sim.system.ny+1:end,imhe));
%             fprintf('Filling the MHE estimation window with data ...\n')
        elseif imhe == sim.config.Ne+1
%             fprintf('Solving for first time the MHE problem ...\n')
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(end-sim.system.ny+1:end,imhe));            
            solve(sim.mhe.nlmhe);                
            sim.data.xestmhe(:,1:sim.config.Ne+1) = sim.mhe.nlmhe.Xtraj;
        else
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(end-sim.system.ny+1:end,imhe));
            solve(sim.mhe.nlmhe);    
            sim.data.xestmhe(:,imhe)     = sim.mhe.nlmhe.x_k;
            %
            sim.data.Jmhe(1,imhe)        = sim.mhe.nlmhe.Jnum;            
            sim.data.Pinv_mhe(:,:,imhe)  = full(inv(sim.mhe.nlmhe.P));
%             compute_remaining_time(sim,imhe,'nlmhe');
        end        
        updateInput(sim.mhe.nlmhe,sim.data.Uk1(1,imhe));
        sim.data.timemhe(1,imhe)     = toc;
        %
        sim.data.max_eig_P_mhe(imhe) = max(real(eig(inv(full(sim.mhe.nlmhe.P)))));
        sim.data.min_eig_P_mhe(imhe) = min(real(eig(inv(full(sim.mhe.nlmhe.P)))));
    end
    % fprintf('\n\n')
    for imheAvg=1:sim.config.N    
        tic;
        yAvg = [];
            for i=1:sim.system.ny
                yAvg = [yAvg;sum(sim.data.ysim1(i:sim.system.ny:end,imheAvg))];
            end
            yAvg = yAvg./sim.config.numSensors;
            %
        if imheAvg < sim.config.Ne+1            
            updateMeasurement(sim.mheAvg.nlmheAvg,yAvg);
%             fprintf('Filling the MHEAVG estimation window with data ...\n')
        elseif imheAvg == sim.config.Ne+1
%             fprintf('Solving for first time the MHEAVG problem ...\n')
            updateMeasurement(sim.mheAvg.nlmheAvg,yAvg);            
            solve(sim.mheAvg.nlmheAvg);            
            sim.data.xestmheAvg(:,1:sim.config.Ne+1) = sim.mheAvg.nlmheAvg.Xtraj;
        else
            updateMeasurement(sim.mheAvg.nlmheAvg,yAvg);
            solve(sim.mheAvg.nlmheAvg);
            sim.data.xestmheAvg(:,imheAvg)     = sim.mheAvg.nlmheAvg.x_k;
            %
            sim.data.JmheAvg(1,imheAvg)        = sim.mheAvg.nlmheAvg.Jnum;            
            sim.data.Pinv_mheAvg(:,:,imheAvg)  = full(inv(sim.mheAvg.nlmheAvg.P));
%             compute_remaining_time(sim,imheAvg,'nlmheavg');
        end        
        updateInput(sim.mheAvg.nlmheAvg,sim.data.Uk1(1,imheAvg));
        sim.data.timemheAvg(1,imheAvg)     = toc;
        %
        sim.data.max_eig_P_mheAvg(imheAvg) = max(real(eig(inv(full(sim.mheAvg.nlmheAvg.P)))));
        sim.data.min_eig_P_mhe(imheAvg) = min(real(eig(inv(full(sim.mheAvg.nlmheAvg.P)))));
    end
    fprintf('\n\n')
    for imheMulY=1:sim.config.N
        tic;
        if imheMulY < sim.config.Ne+1
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
            fprintf('Filling the MHEMULY estimation window with data ...\n')
        elseif imheMulY == sim.config.Ne+1
            fprintf('Solving for first time the MHEMULY problem ...\n')
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
            solve(sim.mheMulY.nlmheMulY);    
            sim.data.xestmheMulY(:,1:sim.config.Ne+1) = sim.mheMulY.nlmheMulY.Xtraj;
        else
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
            solve(sim.mheMulY.nlmheMulY);
            sim.data.xestmheMulY(:,imheMulY) = sim.mheMulY.nlmheMulY.x_k;
            %
            sim.data.JmheMulY(1,imheMulY)    = sim.mheMulY.nlmheMulY.Jnum;
            sim.data.timemheMulY(1,imheMulY) = toc;
            sim.data.Pinv_mheMulY(:,:,imheMulY) = full(inv(sim.mheMulY.nlmheMulY.P));
            compute_remaining_time(sim,imheMulY,'nlmhemuly');
        end        
        updateInput(sim.mheMulY.nlmheMulY,sim.data.Uk1(1,imheMulY));
    end
    fprintf('\n\n')
    for imheConvexY=1:sim.config.N
        tic;
        yaux = [];
        for i=1:sim.system.ny
            yaux = blkdiag(yaux,sim.data.ysim1(i:sim.system.ny:sim.system.ny*sim.config.numSensors,imheConvexY)');
        end
        if imheConvexY < sim.config.Ne+1                           
%             updateMeasurement(sim.mheConvexY.nlmheConvexY,sim.data.ysim1(:,imheConvexY));
updateMeasurement(sim.mheConvexY.nlmheConvexY,yaux);
%             fprintf('Filling the MHECONVEXY estimation window with data ...\n')
        elseif imheConvexY == sim.config.Ne+1
%             fprintf('Solving for first time the MHECONVEXY problem ...\n')
%             updateMeasurement(sim.mheConvexY.nlmheConvexY,reshape(sim.data.ysim1(:,imheConvexY),sim.system.ny,sim.config.numSensors)); 
updateMeasurement(sim.mheConvexY.nlmheConvexY,yaux);
            solve(sim.mheConvexY.nlmheConvexY);                
            sim.data.xestmheConvexY(:,1:sim.config.Ne+1) = sim.mheConvexY.nlmheConvexY.Xtraj;
        else
%             updateMeasurement(sim.mheConvexY.nlmheConvexY,reshape(sim.data.ysim1(:,imheConvexY),sim.system.ny,sim.config.numSensors));
updateMeasurement(sim.mheConvexY.nlmheConvexY,yaux);
            solve(sim.mheConvexY.nlmheConvexY);    
            sim.data.xestmheConvexY(:,imheConvexY) = sim.mheConvexY.nlmheConvexY.x_k;
            %            
            sim.data.JmheConvexY(1,imheConvexY)    = sim.mheConvexY.nlmheConvexY.Jnum;            
            sim.data.Pinv_mheConvexY(:,:,imheConvexY) = full(inv(sim.mheConvexY.nlmheConvexY.P));
%             compute_remaining_time(sim,imheConvexY,'nlmheconvexy');
        end        
        updateInput(sim.mheConvexY.nlmheConvexY,sim.data.Uk1(1,imheConvexY));
        sim.data.timemheConvexY(1,imheConvexY) = toc;
    end    
    for imheLoss=1:sim.config.N      
        tic;
        if imheLoss < sim.config.Ne+1
            updateMeasurement(sim.mheLoss.nlmheLoss,sim.data.ysim1(end-sim.system.ny+1:end,imheLoss));
            % fprintf('Filling the MHELOSS estimation window with data ...\n')
        elseif imheLoss == sim.config.Ne+1
            % fprintf('Solving for first time the MHELOSS problem ...\n')
            updateMeasurement(sim.mheLoss.nlmheLoss,sim.data.ysim1(end-sim.system.ny+1:end,imheLoss));            
            solve(sim.mheLoss.nlmheLoss);                
            sim.data.xestmheLoss(:,1:sim.config.Ne+1) = sim.mheLoss.nlmheLoss.Qtraj;
            sim.data.mheLoss.c = [sim.data.mheLoss.c, sim.mheLoss.nlmheLoss.cvVals];
        else
            updateMeasurement(sim.mheLoss.nlmheLoss,sim.data.ysim1(end-sim.system.ny+1:end,imheLoss));
            solve(sim.mheLoss.nlmheLoss);            
            sim.data.xestmheLoss(:,imheLoss) = sim.mheLoss.nlmheLoss.qk;
            sim.data.mheLoss.c = [sim.data.mheLoss.c, sim.mheLoss.nlmheLoss.cvVals];
            %
%             sim.data.JmheLoss(1,imheLoss) = sim.mheLoss.nlmheLoss.Jnum;            
%             sim.data.Pinv_mheLoss(:,:,imheLoss)  = full(inv(sim.mheLoss.nlmheLoss.P));
            % compute_remaining_time(sim,imheLoss,'nlmheLoss');
        end        
        updateInput(sim.mheLoss.nlmheLoss,sim.data.Uk1(1,imheLoss));
        sim.data.timemheLoss(1,imheLoss) = toc;
        %
%         sim.data.max_eig_P_mheLoss(imheLoss) = max(real(eig(inv(full(sim.mheLoss.nlmheLoss.P)))));
%         sim.data.min_eig_P_mheLoss(imheLoss) = min(real(eig(inv(full(sim.mheLoss.nlmheLoss.P)))));
        %
%         sim.data.mheLoss.alpha_w = [sim.data.mheLoss.alpha_w, sim.mheLoss.nlmheLoss.alpha_w];
%         sim.data.mheLoss.alpha_v = [sim.data.mheLoss.alpha_v, sim.mheLoss.nlmheLoss.alpha_v];
    end
    % fprintf('\n\n')
    % for i=1:sim.config.N
    %     tic;
    %     updateMeasurement(sim.ekf.ekf_joint,sim.data.ysim1(:,i));
    %     solve(sim.ekf.ekf_joint);
    %     sim.data.xestjoint(:,i)     = sim.ekf.ekf_joint.x_k;
    %     sim.data.alphaJoint(:,i)    = sim.ekf.ekf_joint.alpha_k;
    %     sim.data.prmJoint(:,i)      = sim.system.mBA*sim.data.alphaJoint(:,i);
    %     updateInput(sim.ekf.ekf_joint,sim.data.Uk1(1,i));
    %     compute_remaining_time(sim,i,'ekf');
    % end 
    for imheAle=1:sim.config.N    
        tic;
        if imheAle < sim.config.Ne+1
            updateMeasurement(sim.mhe.nlmheAle,sim.data.ysim1(end-sim.system.ny+1:end,imheAle));
%             fprintf('Filling the MHEALE estimation window with data ...\n')
        elseif imheAle == sim.config.Ne+1
%             fprintf('Solving for first time the MHEALE problem ...\n')
            updateMeasurement(sim.mhe.nlmheAle,sim.data.ysim1(end-sim.system.ny+1:end,imheAle));            
            solve(sim.mhe.nlmheAle);                
            sim.data.xestmheAle(:,1:sim.config.Ne+1) = sim.mhe.nlmheAle.Xtraj;
        else
            updateMeasurement(sim.mhe.nlmheAle,sim.data.ysim1(end-sim.system.ny+1:end,imheAle));
            solve(sim.mhe.nlmheAle);    
            sim.data.xestmheAle(:,imheAle)  = sim.mhe.nlmheAle.x_k;
            %            
%             compute_remaining_time(sim,imheAle,'nlmheAle');
        end        
        updateInput(sim.mhe.nlmheAle,sim.data.Uk1(1,imheAle));
        sim.data.timemheAle(1,imheAle)     = toc;
    end
end

function compute_remaining_time(sim,i,mthd)
    time_elapsed        = toc;
    sim.time.total_sim  = sim.time.total_sim + time_elapsed;
    remaining_time      = time_elapsed*(sim.config.N-i);
    remaining_time_min  = floor(remaining_time/60);
    remaining_time_sec  = floor(60*(remaining_time/60 - remaining_time_min));
    %
    fprintf(['Processing ',mthd,' method. ','Remaining time: ', num2str(remaining_time_min), ' [min.], ', num2str(remaining_time_sec), ' [Sec.]\n'])    
end

function sim = compute_errors(sim)
    Ne_plt  = 0;
    t       = linspace((Ne_plt+1)*sim.config.Ts,sim.config.tf,sim.config.N-Ne_plt);
    
%     figure; hold on;
    for k=1:sim.config.numTrials
        norm_mhe        = zeros(1,sim.config.N);
        norm_mheAvg     = zeros(1,sim.config.N);
        norm_mheMulY    = zeros(1,sim.config.N);
        norm_mheConvexY = zeros(1,sim.config.N);
        norm_ekf        = zeros(1,sim.config.N);
        norm_mheLoss    = zeros(1,sim.config.N);
        norm_mheAle     = zeros(1,sim.config.N);
        for i=1:sim.config.N
            norm_mhe(i)     = norm(sim.data.xestmhe_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            norm_mheAvg(i)  = norm(sim.data.xestmheAvg_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            % norm_mheMulY(i)  = norm(sim.data.xestmheMulY_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            norm_mheConvexY(i)  = norm(sim.data.xestmheConvexY_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            norm_mheLoss(i)  = norm(sim.data.xestmheLoss_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            % norm_ekf(i)     = norm(sim.data.xestekf_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            norm_mheAle(i)     = norm(sim.data.xestmheAle_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
        end
        sim.mhe.estErr = mean(norm_mhe);
        sim.mheAvg.estErr = mean(norm_mheAvg);
        sim.mheConvexY.estErr = mean(norm_mheConvexY);
        sim.mheLoss.estErr = mean(norm_mheLoss);
        sim.mheAle.estErr = mean(norm_mheAle);
        %
        sim.data.estErrs = [sim.data.estErrs, [sim.mhe.estErr; sim.mheAvg.estErr; sim.mheConvexY.estErr; sim.mheLoss.estErr; sim.mheAle.estErr]];
%         plot(t(1:end-1),sim.bounds.estimation_error_x_opt+sim.bounds.estimation_error_f_opt,'k','LineWidth',1.5);
%         plot(t,norm_mhe,'r',t,norm_mheAvg,'c',t,norm_mheConvexY,'m',t,norm_mheLoss,'LineWidth',2);grid on;
%         plot(t,norm_mhe,'r',t,norm_mheAvg,'c',t,norm_mheConvexY,'m',t,norm_mheLoss,'y',t,norm_mheAle,'b','LineWidth',2);grid on;
    end
    sim.data.estErrsGlb = [sim.data.estErrsGlb, mean(sim.data.estErrs,2)];
%     grid on;
%     legend({'$MHE$','$MHE_{AVG}$','$MHE_{CONVEXY}$','$MHE_{LOSS}$'},'Interpreter','latex','FontSize',17);
%     legend({'$MHE$','$MHE_{AVG}$','$MHE_{CONVEXY}$','$MHE_{LOSS}$','$MHE_{BAYESIAN}$'},'Interpreter','latex','FontSize',17);
    %
    fprintf('\n');
    fprintf(['err_mhe: ',num2str(sim.data.estErrsGlb(1)),'\n']);
    fprintf(['err_avg: ',num2str(sim.data.estErrsGlb(2)),'\n']);
    fprintf(['err_cvx: ',num2str(sim.data.estErrsGlb(3)),'\n']);
    fprintf(['err_loss: ',num2str(sim.data.estErrsGlb(4)),'\n']);
    fprintf(['err_bayesian: ',num2str(sim.data.estErrsGlb(5)),'\n']);
end

function sim = plt_estimation(sim)
    Ne_plt  = 0;
    t       = linspace((Ne_plt+1)*sim.config.Ts,sim.config.tf,sim.config.N-Ne_plt);
    
    figure; hold on; grid on;    
    for k=1:sim.config.numTrials        
%         plot(t,sim.data.xestmhe_plt((k-1)*sim.system.nx+1,:),'ro',t,sim.data.xestmheAvg_plt((k-1)*sim.system.nx+1,:),'c^',...
%              t,sim.data.xestmheConvexY_plt((k-1)*sim.system.nx+1,:),'md',t,sim.data.xestmheLoss_plt((k-1)*sim.system.nx+1,:),'bv','LineWidth',2);
        plot(t,sim.data.xestmhe_plt((k-1)*sim.system.nx+1,:),'yo',t,sim.data.xestmheAvg_plt((k-1)*sim.system.nx+1,:),'c^',...
             t,sim.data.xestmheConvexY_plt((k-1)*sim.system.nx+1,:),'md',t,sim.data.xestmheLoss_plt((k-1)*sim.system.nx+1,:),'r',...
             t,sim.data.xestmheAle_plt((k-1)*sim.system.nx+1,:),'b','LineWidth',2);
    end
    plot(t,sim.data.xsim1(1,:),'k','LineWidth',2)
    plot(t,sim.data.ysim1(1,:),'g','LineWidth',1)
%     legend({'$MHE$','$MHE_{AVG}$','$MHE_{CONVEXY}$','$MHE_{LOSS}$','$x_1$'},'Interpreter','latex','FontSize',17);
    legend({'$MHE$','$MHE_{AVG}$','$MHE_{CONVEXY}$','$MHE_{LOSS}$','$MHE_{BAYESIAN}$','$x_1$'},'Interpreter','latex','FontSize',17);
    %
    figure; hold on; grid on;    
    for k=1:sim.config.numTrials        
%         plot(t,sim.data.xestmhe_plt(k*sim.system.nx,:),'ro',t,sim.data.xestmheAvg_plt(k*sim.system.nx,:),'c^',...
%              t,sim.data.xestmheConvexY_plt(k*sim.system.nx,:),'md',t,sim.data.xestmheLoss_plt(k*sim.system.nx,:),'bv','LineWidth',2);
        plot(t,sim.data.xestmhe_plt(k*sim.system.nx,:),'yo',t,sim.data.xestmheAvg_plt(k*sim.system.nx,:),'c^',...
             t,sim.data.xestmheConvexY_plt(k*sim.system.nx,:),'md',t,sim.data.xestmheLoss_plt(k*sim.system.nx,:),'r',...
             t,sim.data.xestmheAle_plt(k*sim.system.nx,:),'b','LineWidth',2);
    end
    plot(t,sim.data.xsim1(2,:),'k','LineWidth',2)
    plot(t,sim.data.ysim1(2,:),'g','LineWidth',1)
%     legend({'$MHE$','$MHE_{AVG}$','$MHE_{CONVEXY}$','$MHE_{LOSS}$','$x_2$'},'Interpreter','latex','FontSize',17);
    legend({'$MHE$','$MHE_{AVG}$','$MHE_{CONVEXY}$','$MHE_{LOSS}$','$MHE_{BAYESIAN}$','$x_2$'},'Interpreter','latex','FontSize',17);
end

function sim = plt(sim)
%     sim = compute_errors(sim);
    sim = plt_estimation(sim);
end

function sim = compute_bounds(sim,x0_err,plt)    
    % Compute MHE error bounds ************************************
    iioss_f_bounds = [];
    
    iioss_f_bounds.c_beta_lb = 0;
    iioss_f_bounds.c_beta_ub = 2;
    
    iioss_f_bounds.exp_beta_lb = 0;
    iioss_f_bounds.exp_beta_ub = inf;
    
    iioss_f_bounds.q_lb = 0;
    iioss_f_bounds.q_ub = inf;
    
    iioss_f_bounds.c_gamma_3_lb = 0;
    iioss_f_bounds.c_gamma_3_ub = 1;
    
    iioss_f_bounds.c_gamma_4_lb = 0;
    iioss_f_bounds.c_gamma_4_ub = 1;
    
    iioss_f_bounds.gamma_3_lb = 1;
    iioss_f_bounds.gamma_3_ub = 2;
    
    iioss_f_bounds.gamma_4_lb = 1;
    iioss_f_bounds.gamma_4_ub = 2;
    
    iioss_f_bounds.e_gamma_3_lb = 0;
    iioss_f_bounds.e_gamma_3_ub = inf;
    iioss_f_bounds.e_gamma_4_lb = 0;
    iioss_f_bounds.e_gamma_4_ub = inf;

%             con = [con, exp_beta/2 <= q <= inf];
    %
    open_bounds = [];
    open_bounds.system = 'scalar_1';
    open_bounds.load_iioss = true;
    open_bounds.load_mhe = false;
    open_bounds.load_mhe_opt = true;
    %
            
%     sim.bounds = compute_mhe_bounds(sim.system.nx,sim.system.nw,sim.system.nu,sim.system.ny,sim.system.nv,sim.mhe.mm.Ne,...
%         sim.config.Ts,sim.config.N,inv(full(sim.mhe.mm.P)),inv(full(sim.mhe.mm.Pa)),full(sim.mhe.mm.Qx),full(sim.mhe.mm.Rx),...
%         full(sim.mhe.mm.Ra),full(sim.mhe.mm.Qa),[],sim.estimator.constraints.Xe,[-3.*sim.data.processDistVariance,3.*sim.data.processDistVariance],...
%         sim.estimator.constraints.De,[-3.*sim.data.measNoiseVariance,3.*sim.data.measNoiseVariance],sim.estimator.constraints.Uc,sim.estimator.maxDF,sim.dynamic_anon.f_anon,sim.dynamic.anon.h_anon,x0_err,0,plt);    
    sim.bounds = compute_mhe_bounds(sim.system.nx,sim.system.nw,sim.system.nu,sim.system.ny,sim.system.nv,20,...
        sim.config.Ts,2*sim.config.N,diag(ones(1,sim.system.nx)),diag(ones(1,sim.mhe.mm.Nm)),diag(ones(1,sim.system.nx)),diag(ones(1,sim.system.nv)),...
        diag(ones(1,sim.system.nv)),diag(ones(1,sim.system.nd)),[],sim.estimator.constraints.Xe,[-3.*sim.data.processDistVariance,3.*sim.data.processDistVariance],...
        sim.estimator.constraints.De,[-3.*sim.data.measNoiseVariance,3.*sim.data.measNoiseVariance],sim.estimator.constraints.Uc,sim.estimator.maxDF,sim.dynamic_anon.f_anon,sim.dynamic.anon.h_anon,x0_err,iioss_f_bounds,0,plt,open_bounds);
end
