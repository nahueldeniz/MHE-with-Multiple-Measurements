% #########################################################################
%
% Author: Nestor Nahuel Deniz.
% 
%
% #########################################################################


function sim = husky_mulY()
    clear all;% clc;% close all;
    import casadi.*
    sim = [];
    sim = init(sim);
    sim = gen_constraints(sim);
    sim = gen_state_trajectory(sim);
    sim = gen_x0bar(sim);
    sim = build_estimators(sim);
    sim = set_init_condition_estimators(sim);
    sim = estimate(sim);
    sim = append_to_plot(sim);
    sim = plt(sim);
end

function sim = build_setup(sim)
    fileId = '/home/nahuel/Dropbox/PosDoc/AC3E/MHEMULY/data/cel9-muly.mat';
    sim.filedExpData = load(fileId,'-mat');
    % SIMULATION PARAMETERS ===================================================
    sim.config.t0       = 0;                                                 % initial time of sym. [Seg]
    sim.config.tf       = sum(sim.filedExpData.deltaT);                                                % final time of sym. [Seg]
    sim.config.Ts       = mean(sim.filedExpData.deltaT);                                               % sampling period, 0=>discrete-time system
    %
    sim.mheMulY.numOutputSensors = 3;
    sim.mheConvexY.numOutputSensors = 3;
    % Terraza AC3E
    sim.config.LAT0    = -33.034115;
    sim.config.LON0    = -71.592205;
end

function sim = gen_dynamic(sim)
    % System Parameters _______________________________________________________
    sim.system.nq       = 3;                            % vector state dimension
    sim.system.nw       = 3;                            % process noise input dimension
    sim.system.ny       = 3;                            % output dimension
    sim.system.nu       = 2;                            % input dimension
    sim.system.nv       = 3;
    % Casadi variabes ----------------------------------------------------
    sim.dynamic.q       = casadi.MX.sym('q',sim.system.nq);
    sim.dynamic.w       = casadi.MX.sym('w',sim.system.nw);
    sim.dynamic.u       = casadi.MX.sym('u',sim.system.nu);
    % Dynamic of the system ----------------------------------------------
    sim.dynamic.f_rhs   = [sim.dynamic.u(2) * cos(sim.dynamic.q(3));...
                           sim.dynamic.u(2) * sin(sim.dynamic.q(3));...
                           sim.dynamic.u(1) ];
    sim.dynamic.f       = casadi.Function('f', {sim.dynamic.q,sim.dynamic.w,sim.dynamic.u}, {sim.dynamic.f_rhs});  
    % Output of the system ------------------------------------------------
    sim.dynamic.h_rhs   = sim.dynamic.q;
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
end

function sim = append_to_plot(sim)
    sim.data.xsim_plt       = [sim.data.xsim_plt;sim.data.xsim1];
    sim.data.xestmhe_plt    = [sim.data.xestmhe_plt;sim.data.xestmhe];
    sim.data.xestmheAvg_plt    = [sim.data.xestmheAvg_plt;sim.data.xestmheAvg];
    sim.data.xestmheMulY_plt = [sim.data.xestmheMulY_plt;sim.data.xestmheMulY];
    sim.data.xestmheConvexY_plt = [sim.data.xestmheConvexY_plt;sim.data.xestmheConvexY];
    sim.data.xestekf_plt    = [sim.data.xestekf_plt;sim.data.xestjoint];
    sim.data.xestmheLoss_plt    = [sim.data.xestmheLoss_plt;sim.data.xestmheLoss];
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
end

function sim = reserve_memory(sim)
    % Variables: states, outputs, etc. ____________________________________
    sim.data.xsim1          = [];%zeros(sim.system.nq,sim.config.N+1);
    sim.data.ysim1          = [];%&zeros(sim.system.ny*sim.mheMulY.numOutputSensors,sim.config.N);
    %
    sim.data.xsim_plt       = [];
    sim.data.xestmhe_plt    = [];
    sim.data.xestmheAvg_plt    = [];
    sim.data.xestmheMulY_plt = [];
    sim.data.xestmheConvexY_plt = [];
    sim.data.xestekf_plt    = [];
    sim.data.xestmheLoss_plt    = [];
    %
    sim.data.mheLoss.alpha_w = [];
    sim.data.mheLoss.alpha_v = [];
    %
    sim.data.xestmhe         = [];%NaN(sim.system.nq,sim.config.N);
    sim.data.xestmheAvg      = [];%&NaN(sim.system.nq,sim.config.N);
    sim.data.xestmheMulY     = [];%NaN(sim.system.nq,sim.config.N);
    sim.data.xestmheConvexY  = [];%NaN(sim.system.nq,sim.config.N);
    sim.data.xestjoint       = [];%NaN(sim.system.nq,sim.config.N);
    sim.data.xestmheLoss     = [];%NaN(sim.system.nq,sim.config.N);
    %       
end

function sim = gen_constraints(sim)
    % Constraints _________________________________________________
    sim.estimator.constraints.Xe      = [];
    sim.estimator.constraints.We      = [];
    sim.estimator.constraints.Ve      = [];
    sim.estimator.constraints.Uc      = [];
    % only for the muliple model mhe    
end

function sim = gen_x0bar(sim)
    sim.estimator.x0bar = sim.data.xsim1(1:sim.system.nq,1);
end

function sim = build_nlmhe(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mhe.Ne      = 2;
    sim.mhe.l_x     = 1e3;    
    sim.mhe.Qx      = diag([10 10 1]);                           % states stage-cost
    sim.mhe.P       = sim.mhe.Qx;%sim.mhe.l_x.*casadi.DM.eye(sim.system.nq);        % init. val. of the arrival-cost weight matrix
    sim.mhe.Rx      = diag([15 15 5]);                 % measurement noise stage-cost
    sim.mhe.RuP     = [];%diag([10 200]);                          % measurement noise stage-cost
    sim.mhe.Tsmhe   = sim.config.Ts;                                    % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mhe.nlmhe   = nlmheCasadi(sim.mhe.Ne,sim.dynamic.q,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mhe.Qx,sim.mhe.Rx,sim.mhe.RuP,sim.system.nq,sim.system.nu,sim.system.ny,sim.system.nw,sim.mhe.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'no_adap');
    setSigma(sim.mhe.nlmhe,10*max(0.05,1e-3));
    setC(sim.mhe.nlmhe,2*sim.mhe.l_x);
end

function sim = build_nlmheAvg(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheAvg.Ne      = 2;
    sim.mheAvg.l_x     = 1e3;    
    sim.mheAvg.Qx      = diag([10 10 1]);
    sim.mheAvg.P       = sim.mheAvg.Qx;%sim.mheAvg.l_x.*casadi.DM.eye(sim.system.nq);  % init. val. of the arrival-cost weight matrix
    sim.mheAvg.Rx      = diag([15 15 5]);              % measurement noise stage-cost
    sim.mheAvg.RuP     = [];%1e1.*eye(sim.system.nu);                       % measurement noise stage-cost
    sim.mheAvg.Tsmhe   = sim.config.Ts;                                 % If 0, mhe does not integrate the systems
    % *************************************************************
    sim.mheAvg.nlmheAvg   = nlmheCasadi(sim.mheAvg.Ne,sim.dynamic.q,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheAvg.Qx,sim.mheAvg.Rx,sim.mheAvg.RuP,sim.system.nq,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheAvg.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'no_adap');
    setSigma(sim.mheAvg.nlmheAvg,10*max(0.05,1e-3));
    setC(sim.mheAvg.nlmheAvg,2*sim.mheAvg.l_x);
end

function sim = build_nlmheMulY(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheMulY.Ne      = 2;
    sim.mheMulY.l_x     = 1e3;
    sim.mheMulY.Qx      = diag([10 10 0.1]);
    sim.mheMulY.P       = sim.mhe.Qx;
    sim.mheMulY.Rx      = diag([5 5 5/3]);
    sim.mheMulY.RuP     = [];
    sim.mheMulY.Tsmhe   = sim.config.Ts;
    % *************************************************************
    sim.mheMulY.nlmheMulY = nlmheMulY(sim.mheMulY.Ne,sim.mheMulY.numOutputSensors,sim.dynamic.q,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheMulY.Qx,sim.mheMulY.Rx,sim.mheMulY.RuP,sim.system.nq,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheMulY.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'no_adap');
    setSigma(sim.mheMulY.nlmheMulY,10*max(0.05,1e-3));
    setC(sim.mheMulY.nlmheMulY,2*sim.mheMulY.l_x);
end

function sim = build_nlmheConvexY(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheConvexY.Ne      = 2;
    sim.mheConvexY.l_x     = 1e3;    
    sim.mheConvexY.Qx      = diag([10 10 0.2]);
    sim.mheConvexY.P       = sim.mheConvexY.Qx;
    sim.mheConvexY.Rx      = diag([500 500 110/3]);
    sim.mheConvexY.RuP     = [];
    sim.mheConvexY.Tsmhe   = sim.config.Ts;
    % *************************************************************
    sim.mheConvexY.nlmheConvexY = nlmheConvexY(sim.mheConvexY.Ne,sim.mheConvexY.numOutputSensors,sim.dynamic.q,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheConvexY.Qx,sim.mheConvexY.Rx,sim.mheConvexY.RuP,sim.system.nq,sim.system.nu,sim.system.ny,sim.system.nw,sim.mheConvexY.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'no_adap');
    setSigma(sim.mheConvexY.nlmheConvexY,10*max(0.05,1e-3));
    setC(sim.mheConvexY.nlmheConvexY,2*sim.mheConvexY.l_x);
end

function sim = build_nlmheLoss(sim)
    % ======================= ESTIMATOR ===========================         
    sim.mheLoss.Ne      = 2;
    sim.mheLoss.l_x     = 1e0;
    sim.mheLoss.P       = sim.mhe.l_x.*casadi.DM.eye(sim.system.nq);        % init. val. of the arrival-cost weight matrix
    sim.mheLoss.cw      = 0.01;                                                % Deviation of process disturbances
    sim.mheLoss.Kw      = 0.1;                                               % Weight of alpha_w in cost function
    sim.mheLoss.cv      = 0.15;                                                % Deviation of measurement disturbances
    sim.mheLoss.Kv      = 0.01;                                               % Weight of alpha_v in cost function
    sim.mheLoss.RuP     = 1e1.*eye(sim.system.nu);                          % measurement noise stage-cost
    sim.mheLoss.Tsmhe   = sim.config.Ts;                                    % If 0, mhe does not integrate the systems
    % *************************************************************
% nlmheCasadi_lossFunctions(N,x,w,u,f_rhs,h_rhs,cw,Kalphaw,cv,Kalphav,Ru,nx,nu,ny,nw,P,XLUBounds,WLUBounds,VLUBounds,ULUBounds,x0bar,integ,arrival_cost,arg_opt)    
    sim.mheLoss.nlmheLoss   = nlmheCasadi_lossFunctions(sim.mheLoss.Ne,sim.dynamic.q,sim.dynamic.w,sim.dynamic.u,sim.dynamic.f_rhs,sim.dynamic.h_rhs,...
        sim.mheLoss.cw,sim.mheLoss.Kw,sim.mheLoss.cv,sim.mheLoss.Kv,sim.mhe.RuP,sim.system.nq,sim.system.nu,sim.system.ny,sim.system.nw,sim.mhe.P,sim.estimator.constraints.Xe,...
        sim.estimator.constraints.We,sim.estimator.constraints.Ve,sim.estimator.constraints.Uc,[],sim.config.Ts,'no_adap');
    setSigma(sim.mheLoss.nlmheLoss,10*max(0.05,1e-3));
    setC(sim.mheLoss.nlmheLoss,2*sim.mheLoss.l_x);
end

function sim = build_ekf_joint(sim)
%     sim.ekf.Qekf        = blkdiag(zeros(sim.system.nq),100.*eye(sim.mhe.mm.Nm));
%     sim.ekf.Rekf        = 0.01;
%     sim.ekf.Phat        = 1e-6.*blkdiag(eye(sim.system.nq),eye(sim.mhe.mm.Nm));
%     sim.ekf.ekf_joint   = kalman_x_alpha(sim.mhe.mm.Nm, sim.system.nq,sim.system.nu,sim.system.ny,sim.system.nw,sim.system.nv,...
%         sim.ekf.Phat,sim.mhe.mm.MMAdt,sim.mhe.mm.MMBdt,sim.mhe.mm.MMCdt,sim.mhe.mm.MMaidt,sim.estimator.x0bar,sim.polytope.alpha0,sim.ekf.Qekf,sim.ekf.Rekf);
end

function sim = gen_state_trajectory(sim)        
    rotMtx  = eye(2);   
    xRtk    = [];
    yRtk    = [];
    for i=1:min(length(sim.filedExpData.groundTruthSensors.RTK.Latitude), length(sim.filedExpData.mobileSensors.Orientation))
        [xx, yy]    = latlon2xy(sim.filedExpData.groundTruthSensors.RTK.Latitude(i), sim.filedExpData.groundTruthSensors.RTK.Longitude(i), sim.config.LAT0, sim.config.LON0);
        xy          = rotMtx * [xx*1000; yy*1000];
        xRtk        = [xRtk, xy(1)];
        yRtk        = [yRtk, xy(2)];
    end    
    oriMob          = unwrap(sim.filedExpData.mobileSensors.Orientation(1,:)*pi/180);
    oriMob          = oriMob-oriMob(1);
    %
    oriVec1 = [];
    for i=1:length(sim.filedExpData.groundTruthSensors.VEC1.Orientation)
        rawVec1 = quat2eul([sim.filedExpData.groundTruthSensors.VEC1.Orientation(1,i),sim.filedExpData.groundTruthSensors.VEC1.Orientation(2,i),sim.filedExpData.groundTruthSensors.VEC1.Orientation(3,i),sim.filedExpData.groundTruthSensors.VEC1.Orientation(4,i)]);
        oriVec1 = [oriVec1, rawVec1(3)];
    end
    oriVec1 = unwrap(oriVec1);
    oriVec1 = oriVec1-oriVec1(1);
    %
    xsimAux  = [xRtk;yRtk;oriMob];
    %
    newLength       = floor(length(xsimAux)/sim.mheMulY.numOutputSensors);    
    %
    sim.data.ysim1  = reshape(xsimAux(:,1:newLength*sim.mheMulY.numOutputSensors),sim.mheMulY.numOutputSensors*sim.system.nq,newLength);
    % Here I should choose which sample
    sim.data.xsim1  = [xRtk;yRtk;oriVec1];
    sim.data.xsim1  = downsample(sim.data.xsim1(:,1:newLength*sim.mheMulY.numOutputSensors)',sim.mheMulY.numOutputSensors)';
    %
    linVelVec1 = [];
    for i=1:length(sim.filedExpData.groundTruthSensors.VEC1.nedVel)
        linVelVec1 = [linVelVec1, double(sqrt(sim.filedExpData.groundTruthSensors.VEC1.nedVel(1,i)^2+sim.filedExpData.groundTruthSensors.VEC1.nedVel(2,i)^2))];
    end
    angVelVec1 = double(sim.filedExpData.groundTruthSensors.VEC1.AngularVelocity(3,:));
    sim.data.Uk1 = [angVelVec1; linVelVec1];
    sim.data.Uk1 = full(downsample(sim.data.Uk1(:,1:newLength*sim.mheMulY.numOutputSensors)',sim.mheMulY.numOutputSensors)');
    %
    indx                    = find(isnan(sim.data.xsim1(1,:)));
    sim.data.xsim1(:,indx)  = [];
    sim.data.ysim1(:,indx)  = [];
    sim.data.Uk1(:,indx)    = [];
    %
    sim.config.N = min([length(sim.data.ysim1), length(sim.data.xsim1), length(sim.data.Uk1)]);
end

function sim = estimate(sim)  
    for imhe=1:sim.config.N
        tic;
        if imhe < sim.mhe.Ne+1
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(1:sim.system.ny,imhe));
%             fprintf('Filling the MHE estimation window with data ...\n')
        elseif imhe == sim.mhe.Ne+1
%             fprintf('Solving for first time the MHE problem ...\n')
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(1:sim.system.ny,imhe));
            solve(sim.mhe.nlmhe);    
            sim.data.xestmhe(:,1:sim.mhe.Ne+1) = sim.mhe.nlmhe.Xtraj;
        else
            updateMeasurement(sim.mhe.nlmhe,sim.data.ysim1(1:sim.system.ny,imhe));
            solve(sim.mhe.nlmhe);            
            sim.data.xestmhe(:,imhe)     = sim.mhe.nlmhe.x_k;
        end        
        updateInput(sim.mhe.nlmhe,sim.data.Uk1(:,imhe));
        %
        sim.data.max_eig_P_mhe(imhe) = max(real(eig(inv(full(sim.mhe.nlmhe.P)))));
        sim.data.min_eig_P_mhe(imhe) = min(real(eig(inv(full(sim.mhe.nlmhe.P)))));
    end
%     fprintf('\n')
    for imheAvg=1:sim.config.N
        tic;
        yAvg = [];
            for i=1:sim.system.ny
                yAvg = [yAvg;sum(sim.data.ysim1(i:sim.system.ny:end,imheAvg))];
            end
            yAvg = yAvg./sim.mheMulY.numOutputSensors;
            %
        if imheAvg < sim.mheAvg.Ne+1            
            updateMeasurement(sim.mheAvg.nlmheAvg,yAvg);
%             fprintf('Filling the MHEAVG estimation window with data ...\n')
        elseif imheAvg == sim.mheAvg.Ne+1
%             fprintf('Solving for first time the MHEAVG problem ...\n')
            updateMeasurement(sim.mheAvg.nlmheAvg,yAvg);
            solve(sim.mheAvg.nlmheAvg);
            sim.data.xestmheAvg(:,1:sim.mheAvg.Ne+1) = sim.mheAvg.nlmheAvg.Xtraj;
        else
            updateMeasurement(sim.mheAvg.nlmheAvg,yAvg);
            solve(sim.mheAvg.nlmheAvg);            
            sim.data.xestmheAvg(:,imheAvg)     = sim.mheAvg.nlmheAvg.x_k;
        end        
        updateInput(sim.mheAvg.nlmheAvg,sim.data.Uk1(:,imheAvg));
        %
        sim.data.max_eig_P_mheAvg(imheAvg) = max(real(eig(inv(full(sim.mheAvg.nlmheAvg.P)))));
        sim.data.min_eig_P_mhe(imheAvg) = min(real(eig(inv(full(sim.mheAvg.nlmheAvg.P)))));
    end
%     fprintf('\n')
    for imheMulY=1:sim.config.N
        tic;
        if imheMulY < sim.mheMulY.Ne+1
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
%             fprintf('Filling the MHEMULY estimation window with data ...\n')
        elseif imheMulY == sim.mheMulY.Ne+1
%             fprintf('Solving for first time the MHEMULY problem ...\n')
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
            solve(sim.mheMulY.nlmheMulY);    
            sim.data.xestmheMulY(:,1:sim.mheMulY.Ne+1) = sim.mheMulY.nlmheMulY.Xtraj;
        else
            updateMeasurement(sim.mheMulY.nlmheMulY,sim.data.ysim1(:,imheMulY));
            solve(sim.mheMulY.nlmheMulY);
            sim.data.xestmheMulY(:,imheMulY) = sim.mheMulY.nlmheMulY.x_k;
        end        
        updateInput(sim.mheMulY.nlmheMulY,sim.data.Uk1(:,imheMulY));
    end
%     fprintf('\n')
    for imheConvexY=1:sim.config.N
        tic;
        if imheConvexY < sim.mheConvexY.Ne+1
            updateMeasurement(sim.mheConvexY.nlmheConvexY,sim.data.ysim1(:,imheConvexY));
%             fprintf('Filling the MHECONVEXY estimation window with data ...\n')
        elseif imheConvexY == sim.mheConvexY.Ne+1
%             fprintf('Solving for first time the MHECONVEXY problem ...\n')
            updateMeasurement(sim.mheConvexY.nlmheConvexY,sim.data.ysim1(:,imheConvexY));
            solve(sim.mheConvexY.nlmheConvexY);    
            sim.data.xestmheConvexY(:,1:sim.mheConvexY.Ne+1) = sim.mheConvexY.nlmheConvexY.Xtraj;
        else
            updateMeasurement(sim.mheConvexY.nlmheConvexY,sim.data.ysim1(:,imheConvexY));
            solve(sim.mheConvexY.nlmheConvexY);
            sim.data.xestmheConvexY(:,imheConvexY) = sim.mheConvexY.nlmheConvexY.x_k;
        end        
        updateInput(sim.mheConvexY.nlmheConvexY,sim.data.Uk1(:,imheConvexY));
    end    
%     for imheLoss=1:sim.config.N
%         tic;
%         if imheLoss < sim.mheLoss.Ne+1
%             updateMeasurement(sim.mheLoss.nlmheLoss,sim.data.ysim1(1:sim.system.ny,imheLoss));
%             fprintf('Filling the MHELOSS estimation window with data ...\n')
%         elseif imheLoss == sim.mheLoss.Ne+1
%             fprintf('Solving for first time the MHELOSS problem ...\n')
%             updateMeasurement(sim.mheLoss.nlmheLoss,sim.data.ysim1(1:sim.system.ny,imheLoss));
%             solve(sim.mheLoss.nlmheLoss);    
%             sim.data.xestmheLoss(:,1:sim.mheLoss.Ne+1) = sim.mheLoss.nlmheLoss.Xtraj;
%         else
%             updateMeasurement(sim.mheLoss.nlmheLoss,sim.data.ysim1(1:sim.system.ny,imheLoss));
%             solve(sim.mheLoss.nlmheLoss);            
%             sim.data.xestmheLoss(:,imheLoss) = sim.mheLoss.nlmheLoss.x_k;
%         end        
%         updateInput(sim.mheLoss.nlmheLoss,sim.data.Uk1(:,imheLoss));
%         %
%         sim.data.max_eig_P_mheLoss(imheLoss) = max(real(eig(inv(full(sim.mheLoss.nlmheLoss.P)))));
%         sim.data.min_eig_P_mheLoss(imheLoss) = min(real(eig(inv(full(sim.mheLoss.nlmheLoss.P)))));
%         %
%         sim.data.mheLoss.alpha_w = [sim.data.mheLoss.alpha_w, sim.mheLoss.nlmheLoss.alpha_w];
%         sim.data.mheLoss.alpha_v = [sim.data.mheLoss.alpha_v, sim.mheLoss.nlmheLoss.alpha_v];
%     end
end

function sim = plt_est_error(sim)
    figure;
    subplot(2,1,1); hold on; grid on;
    plot(sim.data.xestmhe(3,:),'y','LineWidth',2);
    plot(sim.data.xestmheAvg(3,:),'c','LineWidth',2);
    plot(sim.data.xestmheConvexY(3,:),'m','LineWidth',2);
    plot(sim.data.xestmheMulY(3,:),'b','LineWidth',2);
    plot(sim.data.xsim1(3,:),'r','LineWidth',2);
    legend('mhe','Avg','Convex','Muly','true')
    %
    for i=1:sim.mheMulY.numOutputSensors
        plot(sim.data.ysim1(3*i,:),'k')
    end
    %
    err_mhe = (sim.data.xestmhe(3,:)-sim.data.xsim1(3,:)).^2;
    err_avg = (sim.data.xestmheAvg(3,:)-sim.data.xsim1(3,:)).^2;
    err_cvx = (sim.data.xestmheConvexY(3,:)-sim.data.xsim1(3,:)).^2;
    err_muly = (sim.data.xestmheMulY(3,:)-sim.data.xsim1(3,:)).^2;
    err_meas = (sim.data.ysim1(3,:)-sim.data.xsim1(3,:)).^2;

    subplot(2,1,2); hold on; grid on;
    plot(err_mhe,'y','LineWidth',2);
    plot(err_avg,'c','LineWidth',2);
    plot(err_cvx,'m','LineWidth',2);
    plot(err_muly,'b','LineWidth',2);
    plot(err_meas,'k','LineWidth',2);
    legend({'$err_{mhe}$','$err_{Avg}$','$err_{Convex}$','$err_{Muly}$','$err_{measurements}$'},'Interpreter','latex','FontSize',20)
    ylim([0 0.4])    
    %
    fprintf('\n');
    fprintf(['err_mhe: ',num2str(mean(err_mhe)),'\n']);
    fprintf(['err_avg: ',num2str(mean(err_avg)),'\n']);
    fprintf(['err_cvx: ',num2str(mean(err_cvx)),'\n']);
    fprintf(['err_muly: ',num2str(mean(err_muly)),'\n']);
    fprintf(['err_meas: ',num2str(mean(err_meas)),'\n']);
end

function sim = plt(sim)
%     sim = plt_config(sim);
%     sim = plt_phase_diagram(sim);
    sim = plt_est_error(sim);
%     sim = plt_estimation(sim);
end

