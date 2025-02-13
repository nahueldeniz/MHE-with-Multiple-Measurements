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
    %
    sim = gen_x0(sim);    
    %
%     sim = find_max_init_x0_err(sim);
    sim = gen_constraints(sim);
    %    
%     sim = compute_bounds(sim,sim.config.Max_err_x0,false);%sim.estimator.x0_max_error
%     %
    sim = gen_noise_sequences(sim);
    sim = gen_state_trajectory(sim);
    %    
    for i=1:sim.config.num_sim
        sim = gen_x0bar(sim);
        %                
        sim = build_estimators(sim);
        %
        sim = set_init_condition_estimators(sim);
        %
        sim = estimate(sim);
        %
        sim = append_to_plot(sim);
    end    
    %
%     sim = save_vars(sim);
    %
    sim = plt(sim);
end

function sim = build_setup(sim)
    % SIMULATION PARAMETERS ===================================================
    sim.config.t0       = 0;                                                 % initial time of sym. [Seg]
    sim.config.tf       = 25;                                                % final time of sym. [Seg]
    sim.config.Ts       = 0.1;                                               % sampling period, 0=>discrete-time system
    sim.config.N        = (sim.config.tf-sim.config.t0)/sim.config.Ts-1;     % number of steps of sym.
    sim.config.same_seed = false;
    % Noise amplitudes __________________________________________________________
    sim.data.Sw         = 0.15;                                             % process noise amplitude
    sim.data.Sv         = 0.05;                                             % measurement noise amplitude    
    sim.data.outlier    = false;
    %
    sim.config.Max_err_x0 = 1;                                               % deviation of the initial guess
    %
    sim.time.total_sim  = 0;
    %
    sim.config.num_sim  = 1;
    %
    sim.mheMulY.numOutputSensors = 10;
end

function sim = gen_dynamic(sim)
    % System Parameters _______________________________________________________
    sim.system.nx      = 1;                            % vector state dimension
    sim.system.nw      = 1;                            % process noise input dimension
    sim.system.nd      = 1;                            % model mismatch disturbance
    sim.system.ny      = 1;                            % output dimension
    sim.system.nu      = 1;                            % input dimension
    sim.system.nv      = sim.system.ny;
    % Casadi variabes ----------------------------------------------------
    sim.dynamic.x      = casadi.MX.sym('x',sim.system.nx);
    sim.dynamic.w      = casadi.MX.sym('w',sim.system.nw);
    sim.dynamic.u      = casadi.MX.sym('u',sim.system.nu);
    % System's parameters ------------------------------------------------
    sim.system.p1      = 0.5;
    sim.system.p2      = 2;
    % Dynamic of the system ----------------------------------------------
    sim.dynamic.f_rhs  = -1 * sim.system.p1*sim.dynamic.x^sim.system.p2 + sim.dynamic.w + sim.dynamic.u;
    sim.dynamic.f      = casadi.Function('f', {sim.dynamic.x,sim.dynamic.w,sim.dynamic.u}, {sim.dynamic.f_rhs});
    sim.dynamic.fToLin = sim.dynamic.f;
    % Output of the system ------------------------------------------------
    sim.dynamic.h_rhs  = sim.dynamic.x^(1);
    sim.dynamic.h      = casadi.Function('h', {sim.dynamic.x}, {sim.dynamic.h_rhs});        
    % RK4 -----------------------------------------------------------------
    k1      = sim.dynamic.f(sim.dynamic.x,sim.dynamic.w,sim.dynamic.u);
    k2      = sim.dynamic.f(sim.dynamic.x + sim.config.Ts / 2 * k1, sim.dynamic.w, sim.dynamic.u);
    k3      = sim.dynamic.f(sim.dynamic.x + sim.config.Ts / 2 * k2, sim.dynamic.w, sim.dynamic.u);
    k4      = sim.dynamic.f(sim.dynamic.x + sim.config.Ts * k3, sim.dynamic.w, sim.dynamic.u);
    x_rk4   = sim.dynamic.x + sim.config.Ts / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    sim.dynamic.F = casadi.Function('F', {sim.dynamic.x, sim.dynamic.w, sim.dynamic.u}, {x_rk4});    
    % Used only with yhe bounds calculator --------------------------------
    sim.dynamic_anon.f_anon  = @(t,x,w,u)( -1 * sim.system.p1*x^sim.system.p2 + w + u);
    sim.dynamic.anon.h_anon  = @(x,v)(x^(1));
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
    sim = build_ekf_joint(sim);
end

function sim = append_to_plot(sim)
    sim.data.xsim_plt       = [sim.data.xsim_plt;sim.data.xsim1];
    sim.data.xestmhe_plt    = [sim.data.xestmhe_plt;sim.data.xestmhe];
    sim.data.xestmheAvg_plt    = [sim.data.xestmheAvg_plt;sim.data.xestmheAvg];
    sim.data.xestmheMulY_plt = [sim.data.xestmheMulY_plt;sim.data.xestmheMulY];
    sim.data.xestekf_plt    = [sim.data.xestekf_plt;sim.data.xestjoint];
end

function sim = set_init_condition_estimators(sim)
%     set_x0bar(sim.ekf.ekf_joint, sim.estimator.x0bar);
%     set_alpha0bar(sim.ekf.ekf_joint, sim.polytope.alpha0bar);
    set_x0bar(sim.mhe.nlmhe, sim.estimator.x0bar);
    %
    set_x0bar(sim.mheAvg.nlmheAvg, sim.estimator.x0bar);
    %
    set_x0bar(sim.mheMulY.nlmheMulY, sim.estimator.x0bar);
end

function sim = reserve_memory(sim)
    % Number of simulation to carry-on ____________________________________
    sim.data.Jmhe           = casadi.DM.zeros(1,sim.config.N);
    sim.data.JmheAvg        = casadi.DM.zeros(1,sim.config.N);
    sim.data.JmheMulY       = casadi.DM.zeros(1,sim.config.N);
    sim.data.JMMmhe         = [];
    sim.data.arrCost        = [];
    % Variables: states, outputs, etc. ____________________________________
    sim.data.xsim1          = zeros(sim.system.nx,sim.config.N+1);
    sim.data.ysim1          = zeros(sim.system.ny*sim.mheMulY.numOutputSensors,sim.config.N);
    sim.data.Uk1            = zeros(1,sim.config.N);%0.1*idinput([20,1,(sim.config.N+1)/5])';
    %
    sim.data.xsim_plt       = [];
    sim.data.xestmhe_plt    = [];
    sim.data.xestmheAvg_plt    = [];
    sim.data.xestmheMulY_plt = [];
    sim.data.xestekf_plt    = [];
    %
    sim.data.xestmhe         = NaN(sim.system.nx,sim.config.N);
    sim.data.xestmheAvg      = NaN(sim.system.nx,sim.config.N);
    sim.data.xestmheMulY     = NaN(sim.system.nx,sim.config.N);
    sim.data.xestjoint       = NaN(sim.system.nx,sim.config.N);
    %
    sim.data.timemhe         = NaN(1,sim.config.N);
    sim.data.timemheAvg      = NaN(1,sim.config.N);
    sim.data.timemheMulY     = NaN(1,sim.config.N);
    sim.data.timeJoint       = NaN(1,sim.config.N);
    %
    sim.data.Pinv_mhe        = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    sim.data.Pinv_mheAvg     = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    sim.data.Pinv_mheMulY    = NaN(sim.system.nx,sim.system.nx,sim.config.N);
    %
    sim.data.max_eig_P_mhe   = NaN(1,sim.config.N);
    sim.data.max_eig_P_mheAvg = NaN(1,sim.config.N);
    sim.data.max_eig_P_mheMulY = NaN(1,sim.config.N);
    %
    sim.data.min_eig_P_mhe   = NaN(1,sim.config.N);
    sim.data.min_eig_P_mheAvg = NaN(1,sim.config.N);
    sim.data.min_eig_P_mheMulY = NaN(1,sim.config.N);
    %       
end

function sim = gen_constraints(sim)
    % Constraints _________________________________________________
    sim.estimator.constraints.Xe      = [-5 5];
    sim.estimator.constraints.We      = 3.*[-sim.data.Sw sim.data.Sw];
    sim.estimator.constraints.Ve      = 3.*[-sim.data.Sv sim.data.Sv];
    sim.estimator.constraints.Uc      = [];
    % only for the muliple model mhe    
end

function sim = gen_x0(sim)
    sim.system.x0       = 2.5;
    sim.data.xsim1(:,1) = sim.system.x0;
end

function sim = gen_x0bar(sim)
    sim.estimator.x0bar = (sim.system.x0 + (sim.config.Max_err_x0.*(2.*rand(sim.system.nx,1)-1))./sim.system.nx);
end

function sim = build_nlmhe(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mhe.Ne      = 15;
    sim.mhe.l_x     = 1e4;
    sim.mhe.P       = sim.mhe.l_x.*casadi.DM.eye(sim.system.nx);              % init. val. of the arrival-cost weight matrix
    sim.mhe.Qx      = diag(10);                % states stage-cost
    sim.mhe.Rx      = 15.*casadi.DM.eye(sim.system.ny);            % measurement noise stage-cost
    sim.mhe.RuP     = 1e6.*eye(sim.system.nu);                 % measurement noise stage-cost
    sim.mhe.Tsmhe   = sim.config.Ts;                           % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mhe.nlmhe   = nlmheCasadi(sim.mhe.Ne,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mhe.Qx,sim.mhe.Rx,sim.mhe.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mhe.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
    setSigma(sim.mhe.nlmhe,10*max(max(sim.data.Sv),1e-3));
    setC(sim.mhe.nlmhe,2*sim.mhe.l_x);
end

function sim = build_nlmheAvg(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheAvg.Ne      = 15;
    sim.mheAvg.l_x     = 1e4;
    sim.mheAvg.P       = sim.mheAvg.l_x.*casadi.DM.eye(sim.system.nx);              % init. val. of the arrival-cost weight matrix
    sim.mheAvg.Qx      = diag(10);                % states stage-cost
    sim.mheAvg.Rx      = 15.*casadi.DM.eye(sim.system.ny);            % measurement noise stage-cost
    sim.mheAvg.RuP     = 1e6.*eye(sim.system.nu);                 % measurement noise stage-cost
    sim.mheAvg.Tsmhe   = sim.config.Ts;                           % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mheAvg.nlmheAvg   = nlmheCasadi(sim.mheAvg.Ne,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheAvg.Qx,sim.mheAvg.Rx,sim.mheAvg.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheAvg.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
    setSigma(sim.mheAvg.nlmheAvg,10*max(max(sim.data.Sv),1e-3));
    setC(sim.mheAvg.nlmheAvg,2*sim.mheAvg.l_x);
end

function sim = build_nlmheMulY(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheMulY.Ne      = 15;
    sim.mheMulY.l_x     = 1e4;
    sim.mheMulY.P       = sim.mheMulY.l_x.*casadi.DM.eye(sim.system.nx);              % init. val. of the arrival-cost weight matrix
    sim.mheMulY.Qx      = diag(10);                % states stage-cost
    sim.mheMulY.Rx      = 15.*casadi.DM.eye(sim.system.ny);            % measurement noise stage-cost
    sim.mheMulY.RuP     = 1e6.*eye(sim.system.nu);                 % measurement noise stage-cost
    sim.mheMulY.Tsmhe   = sim.config.Ts;                           % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mheMulY.nlmheMulY = nlmheMulY(sim.mheMulY.Ne,sim.mheMulY.numOutputSensors,sim.dynamic.x,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheMulY.Qx,sim.mheMulY.Rx,sim.mheMulY.RuP,sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheMulY.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'adap');
    setSigma(sim.mheMulY.nlmheMulY,10*max(max(sim.data.Sv),1e-3));
    setC(sim.mheMulY.nlmheMulY,2*sim.mheMulY.l_x);
end

function sim = build_ekf_joint(sim)
%     sim.ekf.Qekf        = blkdiag(zeros(sim.system.nx),100.*eye(sim.mhe.mm.Nm));
%     sim.ekf.Rekf        = 0.01;
%     sim.ekf.Phat        = 1e-6.*blkdiag(eye(sim.system.nx),eye(sim.mhe.mm.Nm));
%     sim.ekf.ekf_joint   = kalman_x_alpha(sim.mhe.mm.Nm, sim.system.nx,sim.system.nu,sim.system.ny,sim.system.nw,sim.system.nv,...
%         sim.ekf.Phat,sim.mhe.mm.MMAdt,sim.mhe.mm.MMBdt,sim.mhe.mm.MMCdt,sim.mhe.mm.MMaidt,sim.estimator.x0bar,sim.polytope.alpha0,sim.ekf.Qekf,sim.ekf.Rekf);
end

function sim = gen_noise_sequences(sim)    
    % Generate sequences of noises ________________________________
    if sim.config.same_seed
        sim.config_seed.s       = rng;                          % seed number
        sim.config.seed.ran     = 564967342;                           % set seed for rand generator
        rng(sim.config.seed.ran, 'twister');
    end
    sim.data.wc      = [];
    sim.data.vc      = [];
    for j=1:sim.system.nw
        sim.data.wc   = [sim.data.wc;sim.data.Sw(j).*randn(1,sim.config.N)];
    end
    for i=1:sim.mheMulY.numOutputSensors
        for j=1:sim.system.nv
            sim.data.vc   = [sim.data.vc;sim.data.Sv(j).*randn(1,sim.config.N)];
            if sim.data.outlier
                Ne = 15;%min([sim.mhe.Ne,sim.mheAvg.Ne,sim.mheMulY.Ne]);
                for k=Ne+1:sim.config.N
                    if min(abs(sim.data.vc(end,k-Ne:k))) < 3*sim.data.Sv
                        sim.data.vc(end,k+1-randi(Ne)) = 5*sim.data.vc(end,k+1-randi(Ne+1));
                    end
                end
            end
        end
    end
end

function sim = gen_state_trajectory(sim)    
%     sim.data.prmEst(:,1)     = sim.system.mBA*sim.mhe.mm.alpha0bar;
%     sim.data.prmJoint(:,1)   = sim.system.mBA*sim.mhe.mm.alpha0bar;
    
    for i=1:sim.config.N
        % States __________________________________________________
        sim.data.xsim1(:,i+1) = full(sim.dynamic.F(sim.data.xsim1(:,i),sim.data.wc(:,i),sim.data.Uk1(1,i))); 
        % Measurements ____________________________________________
        for j=1:sim.mheMulY.numOutputSensors
            sim.data.ysim1((j-1)*sim.system.ny+1:j*sim.system.ny,i) = full(sim.dynamic.h(sim.data.xsim1(:,i))) + sim.data.vc((j-1)*sim.system.ny+1:j*sim.system.ny,i);
        end
    end  
    sim.data.xsim1 = sim.data.xsim1(:,1:sim.config.N);
end

function sim = estimate(sim)  
    for imhe=1:sim.config.N
        tic;
        if imhe < sim.mhe.Ne+1
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(1:sim.system.ny,imhe));
            fprintf('Filling the MHE estimation window with data ...\n')
        elseif imhe == sim.mhe.Ne+1
            fprintf('Solving for first time the MHE problem ...\n')
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(1:sim.system.ny,imhe));
            solve(sim.mhe.nlmhe);    
            sim.data.xestmhe(:,1:sim.mhe.Ne+1) = sim.mhe.nlmhe.Xtraj;
        else
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(1:sim.system.ny,imhe));
            solve(sim.mhe.nlmhe);            
            sim.data.xestmhe(:,imhe)     = sim.mhe.nlmhe.x_k;
            %
            sim.data.Jmhe(1,imhe)        = sim.mhe.nlmhe.Jnum;
            sim.data.timemhe(1,imhe)     = toc;
            sim.data.Pinv_mhe(:,:,imhe)  = full(inv(sim.mhe.nlmhe.P));
            compute_remaining_time(sim,imhe,'nlmhe');
        end        
        updateInput(sim.mhe.nlmhe,sim.data.Uk1(1,imhe));
        %
        sim.data.max_eig_P_mhe(imhe) = max(diag(real(eig(inv(full(sim.mhe.nlmhe.P))))));
        sim.data.min_eig_P_mhe(imhe) = min(diag(real(eig(inv(full(sim.mhe.nlmhe.P))))));
    end
    fprintf('\n\n')
    for imheAvg=1:sim.config.N
        tic;
        if imheAvg < sim.mheAvg.Ne+1
            updateMeasurement(sim.mheAvg.nlmheAvg,sum(sim.data.ysim1(:,imheAvg)/sim.mheMulY.numOutputSensors));
            fprintf('Filling the MHEAVG estimation window with data ...\n')
        elseif imheAvg == sim.mheAvg.Ne+1
            fprintf('Solving for first time the MHEAVG problem ...\n')
            updateMeasurement(sim.mheAvg.nlmheAvg,sum(sim.data.ysim1(:,imheAvg)/sim.mheMulY.numOutputSensors));
            solve(sim.mheAvg.nlmheAvg);
            sim.data.xestmheAvg(:,1:sim.mheAvg.Ne+1) = sim.mheAvg.nlmheAvg.Xtraj;
        else
            updateMeasurement(sim.mheAvg.nlmheAvg,sum(sim.data.ysim1(:,imheAvg)/sim.mheMulY.numOutputSensors));
            solve(sim.mheAvg.nlmheAvg);            
            sim.data.xestmheAvg(:,imheAvg)     = sim.mheAvg.nlmheAvg.x_k;
            %
            sim.data.JmheAvg(1,imheAvg)        = sim.mheAvg.nlmheAvg.Jnum;
            sim.data.timemheAvg(1,imheAvg)     = toc;
            sim.data.Pinv_mheAvg(:,:,imheAvg)  = full(inv(sim.mheAvg.nlmheAvg.P));
            compute_remaining_time(sim,imheAvg,'nlmheavg');
        end        
        updateInput(sim.mheAvg.nlmheAvg,sim.data.Uk1(1,imheAvg));
        %
        sim.data.max_eig_P_mheAvg(imheAvg) = max(diag(real(eig(inv(full(sim.mheAvg.nlmheAvg.P))))));
        sim.data.min_eig_P_mhe(imheAvg) = min(diag(real(eig(inv(full(sim.mheAvg.nlmheAvg.P))))));
    end
    fprintf('\n\n')
    for imheMulY=1:sim.config.N
        tic;
        if imheMulY < sim.mheMulY.Ne+1
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
            fprintf('Filling the MHEMULY estimation window with data ...\n')
        elseif imheMulY == sim.mheMulY.Ne+1
            fprintf('Solving for first time the MHEMULY problem ...\n')
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
            solve(sim.mheMulY.nlmheMulY);    
            sim.data.xestmheMulY(:,1:sim.mheMulY.Ne+1) = sim.mheMulY.nlmheMulY.Xtraj;
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
%     fprintf('\n\n')
%     for i=1:sim.config.N
%         tic;
%         updateMeasurement(sim.ekf.ekf_joint,sim.data.ysim1(:,i));
%         solve(sim.ekf.ekf_joint);
%         sim.data.xestjoint(:,i)     = sim.ekf.ekf_joint.x_k;
%         sim.data.alphaJoint(:,i)    = sim.ekf.ekf_joint.alpha_k;
%         sim.data.prmJoint(:,i)      = sim.system.mBA*sim.data.alphaJoint(:,i);
%         updateInput(sim.ekf.ekf_joint,sim.data.Uk1(1,i));
%         compute_remaining_time(sim,i,'ekf');
%     end 
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

function sim = plt_config(sim)
    sim.plt.plt_x_true  = 1;
    sim.plt.plt_x_mhe   = 1;
    sim.plt.plt_x_mheMulY = 1;
    sim.plt.plt_x_ekf   = 1;
end

function sim = plt_est_error(sim)
    Ne_plt  = 0;
    t       = linspace((Ne_plt+1)*sim.config.Ts,sim.config.tf,sim.config.N-Ne_plt);
%     t = 1:sim.config.N;
    f  = figure('Name','err_x');hold on;
%     f1 = figure('Name','err_x1');hold on;
%     f2 = figure('Name','err_x2');hold on;
    
    norm_x      = zeros(1,sim.config.N);    
    
    figure(f);
    hold on;
    for k=1:sim.config.num_sim
        norm_mhe     = zeros(1,sim.config.N);
        norm_mheAvg  = zeros(1,sim.config.N);
        norm_mheMulY = zeros(1,sim.config.N);
        norm_ekf     = zeros(1,sim.config.N);
        for i=1:sim.config.N;
            norm_mhe(i)     = norm(sim.data.xestmhe_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            norm_mheAvg(i)  = norm(sim.data.xestmheAvg_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
            norm_mheMulY(i)  = norm(sim.data.xestmheMulY_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));       
            norm_ekf(i)     = norm(sim.data.xestekf_plt((k-1)*sim.system.nx+1:k*sim.system.nx,i)-sim.data.xsim1(:,i));
        end
%         plot(t(1:end-1),sim.bounds.estimation_error_x_opt+sim.bounds.estimation_error_f_opt,'k','LineWidth',1.5);
        plot(t,norm_mhe,'r',t,norm_mheAvg,'g',t,norm_mheMulY,'b',t,norm_ekf,'y');
    end
%     ylim([0 5])
    grid on;
end

function sim = plt_estimation(sim)
    Ne_plt  = 0;
    t       = linspace((Ne_plt+1)*sim.config.Ts,sim.config.tf,sim.config.N-Ne_plt);
%     t = 1:sim.config.N;
    f  = figure('Name','estimation');hold on;
%     f1 = figure('Name','err_x1');hold on;
%     f2 = figure('Name','err_x2');hold on;
    
    figure(f);
    hold on;
    plot(t,sim.data.xsim1,'k','LineWidth',2)
    for k=1:sim.config.num_sim        
        plot(t,sim.data.xestmhe(k,:),'r',t,sim.data.xestmheAvg_plt,'g',t,sim.data.xestmheMulY_plt(k,:),'b-.',t,sim.data.xestekf_plt(k,:),'y','LineWidth',2);
    end
%     ylim([0 5])
    grid on;
end

function sim = plt(sim)
    sim = plt_config(sim);
%     sim = plt_phase_diagram(sim);
    sim = plt_est_error(sim);
    sim = plt_estimation(sim);
end

function sim = save_vars(sim)
    save(['Ne=',num2str(sim.mhe.Ne),'-Sw=',num2str(sim.data.Sw(1)),'-Sv=',num2str(sim.data.Sv(1)),'-epsilon=',num2str(sim.system.epsilon),'-NumSim=',num2str(sim.config.num_sim),'_vanderpol.mat']);
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
%         full(sim.mhe.mm.Ra),full(sim.mhe.mm.Qa),[],sim.estimator.constraints.Xe,[-3.*sim.data.Sw,3.*sim.data.Sw],...
%         sim.estimator.constraints.De,[-3.*sim.data.Sv,3.*sim.data.Sv],sim.estimator.constraints.Uc,sim.estimator.maxDF,sim.dynamic_anon.f_anon,sim.dynamic.anon.h_anon,x0_err,0,plt);    
    sim.bounds = compute_mhe_bounds(sim.system.nx,sim.system.nw,sim.system.nu,sim.system.ny,sim.system.nv,20,...
        sim.config.Ts,2*sim.config.N,diag(ones(1,sim.system.nx)),diag(ones(1,sim.mhe.mm.Nm)),diag(ones(1,sim.system.nx)),diag(ones(1,sim.system.nv)),...
        diag(ones(1,sim.system.nv)),diag(ones(1,sim.system.nd)),[],sim.estimator.constraints.Xe,[-3.*sim.data.Sw,3.*sim.data.Sw],...
        sim.estimator.constraints.De,[-3.*sim.data.Sv,3.*sim.data.Sv],sim.estimator.constraints.Uc,sim.estimator.maxDF,sim.dynamic_anon.f_anon,sim.dynamic.anon.h_anon,x0_err,iioss_f_bounds,0,plt,open_bounds);
end
