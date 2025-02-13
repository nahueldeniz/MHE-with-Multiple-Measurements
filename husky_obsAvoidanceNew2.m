% #########################################################################
% 
% Obstacle Avoidance for Generalised N-Trailer Systems
% Author: Nestor. N. Deniz - 2024
%
% #########################################################################notification icon


function S = N_trailersGPSOutliers()
    clear all; clc;
    import casadi.*
    % ---------------------------------------------------------------------
    S = init();    
    % Init ROS's things ---------------------------------------------------
    if S.config.dataOnline == true && ~S.config.stab.compFlg; figure('units','normalized','outerposition',[0 0 1 1]); hold on; end
    if ~S.config.SIM == true; S = ROS(S); end
    for num_sims = 1:S.config.NUM_SIMS
        S.data.num_sims = num_sims;
        S = call_init_functions(S);
        %
        if S.config.stab.compFlg
            S = makeSolverFindLipschitzNcConstant(S);
            S = setNc(S);
            %
            % Gen set of final conditions qNc
            %
            nq      = S.system.nq;
            nu      = S.system.nu;
            Nc      = S.config.Nc;
            Nt      = S.config.Nt;
            %
            xRange      = (S.config.stab.xMin:0.01:S.config.stab.xMax);
            yRange      = (S.config.stab.yMin:0.01:S.config.stab.yMax);
            theta0Range = (S.config.stab.thetaMin:1*pi/180:S.config.stab.thetaMax);
            theta1Range = [-pi/2 -pi/4 0 pi/4 pi/2];
            indx        = 1;
            numPoints   = length(xRange)*length(yRange)*length(theta0Range)*length(theta1Range);
            QNc         = zeros(nq,numPoints);
            for x=1:length(xRange)
                for y=1:length(yRange)
                    for theta0=1:length(theta0Range)
                        for theta1=1:length(theta1Range)
                            S           = gen_x0(S,[xRange(x),yRange(y),theta0Range(theta0),theta0Range(theta0)+theta1Range(theta1)]);
                            QNc(:,indx) = S.init_condition.x0;
                            indx = indx+1;
                        end
                    end
                end        
            end  
    %
%     S.initConditionsStab = gpuArray(S.initConditionsStab);
%             S.config.stab.nx0Stab = indx-1;
            %            
            du_lb   = S.mpc.box_constraints.dUluBounds(:,1);
            du_ub   = S.mpc.box_constraints.dUluBounds(:,2);
            u_lb    = S.mpc.box_constraints.UluBounds(:,1);
            u_ub    = S.mpc.box_constraints.UluBounds(:,2);
            q_lb    = S.mpc.box_constraints.QluBounds(:,1);
            q_ub    = S.mpc.box_constraints.QluBounds(:,2);
            qN_lb   = S.mpc.box_constraints.QNluBounds(:,1);
            qN_ub   = S.mpc.box_constraints.QNluBounds(:,2);
            FNt     = S.dynamic.FNt;

            qNc     = casadi.MX.sym('qNc',nq);
            Qk      = {};
            Uk      = {};
            for i=1:Nc
                Qk  = [Qk(:)', {casadi.MX.sym(['Q_{k+' num2str(i),'}'],nq)}];
                if i<Nc
                    Uk  = [Uk(:)', {casadi.MX.sym(['U_{k+' num2str(i-1),'}'],nu)}];
                end
            end
            % Contraints of the optimisation problem
            stateConstraints                = {};
            state_constraints_lb            = casadi.DM.zeros(Nc*nq);
            state_constraints_ub            = casadi.DM.zeros(Nc*nq);
            %
            deltaUConstraints               = {Uk{1}};
            deltaUConstraints_lb            = [[0;0]; repmat(du_lb,Nc-2,1)];
            deltaUConstraints_ub            = [[0;0]; repmat(du_ub,Nc-2,1)];
            %
            for i=1:Nc-1
                Fk = FNt(Qk{i},Uk{i});
                %                
                stateConstraints = [stateConstraints {Qk{i+1}-Fk}];
                if i>1
                    deltaUConstraints = [deltaUConstraints {Uk{i}-Uk{i-1}}];
                end
            end
            stateConstraints = [stateConstraints {Qk{Nc}-qNc}];
            %
            J               = -(1/2)*mtimes((Qk{1}([Nt+1,2*Nt+2:2*Nt+3])-qNc([Nt+1,2*Nt+2:2*Nt+3]))',(Qk{1}([Nt+1,2*Nt+2:2*Nt+3])-qNc([Nt+1,2*Nt+2:2*Nt+3])));
            %
            optVar          = [Qk(:)',Uk(:)'];
            optVar_lb   = vertcat( repmat(q_lb, Nc-1, 1),...    q_k
                                              qN_lb,...                 q_N
                                       repmat(u_lb, Nc-1, 1));%       u_k
            optVar_ub   = vertcat( repmat(q_ub, Nc-1, 1),...    q_k
                                              qN_ub,...                 q_N
                                       repmat(u_ub, Nc-1, 1));%       u_k
            optParam        = {qNc};
            nlp_constraints         = [stateConstraints(:)', deltaUConstraints(:)'];
            nlp_constraints_lb      = [state_constraints_lb(:); deltaUConstraints_lb(:)];
            nlp_constraints_ub      = [state_constraints_ub(:); deltaUConstraints_ub(:)];
            %
            problem = struct('f',J,'x',vertcat(optVar{:}),'g',vertcat(nlp_constraints{:}),'p',vertcat(optParam{:}));
            nlpoptions                                  = struct;
            nlpoptions.ipopt.max_iter                   = 2000;  %2000
            nlpoptions.ipopt.print_level                = 0;
            nlpoptions.print_time                       = 0;
            nlpoptions.ipopt.acceptable_tol             = 1e-8;
            nlpoptions.ipopt.acceptable_obj_change_tol  = 1e-6;
            solver  = casadi.nlpsol('solver', 'ipopt', problem, nlpoptions);
            %
            trajStab    = zeros(nq,Nc,numPoints);
            UStab       = zeros(nu,Nc-1,numPoints);
            parfor i=1:numPoints
                sol              = solver('x0' ,[],'lbx',optVar_lb,'ubx',optVar_ub,'lbg',nlp_constraints_lb,'ubg',nlp_constraints_ub,'p',QNc(:,i));             
                x_opt            = full(sol.x);
                trajStab(:,:,i)  = reshape(x_opt(1:nq*Nc),nq,Nc);
                UStab(:,:,i)     = reshape(x_opt(nq*Nc+1:end),nu,Nc-1);
            end
            S.trajStab  = trajStab;
            S.UStab     = UStab;

        else
        %
            S.config.iters = 1;
            while (S.config.time(end) < S.config.tf)% && ~(S.path.reach_end_mhempc || ~S.config.mpc)
                if check_end_condition(S); break; else
                    S = CALL_SIMULATOR(S);
                end
                % #############################################################        
                S.config.time   = [S.config.time, S.config.time(end)+S.config.Ts];
                S.config.iters  = S.config.iters+1;
                a               = [num_sims, S.config.iters, S.exec_time.t_tot(end)/S.config.Ts]
                % -------------------------------------------------------------
                if S.config.dataOnline == true
                    tic;
                    plot_mono2(S, S.data.xest(:,end));
                    if S.config.video
                        frame = getframe();
                        writeVideo(S.video.writerObj, frame)
                    end
                    S.exec_time.t_plt = [S.exec_time.t_plt, toc];
                end
                S.exec_time.t_mis = [S.exec_time.t_mis, toc];
                % -------------------------------------------------------------
                if S.config.mpc == true
                    S.data.mhempc.performance.xfut{num_sims, S.config.iters-1}  = S.mpc.mpcCasadi.Qtraj;
                    S.data.mhempc.performance.J{num_sims, S.config.iters-1}     = S.mpc.mpcCasadi.Jnum;
                    S.data.mhempc.performance.qref{num_sims, S.config.iters-1}  = S.mpc.mpcCasadi.Qtraj;
                    S.data.mhempc.performance.XYR{num_sims, S.config.iters-1}   = S.mpc.mpcCasadi.XYRobst;
                end
            end
            if ~S.config.SIM
                S = write_ctrl(S, 0, 0); % Stop the vehiclee            
            end
            if S.config.video
                close(S.video.writerObj);
            end
            % -----------------------------------------------------------------
            S.data.mhempc.performance.xsim{num_sims} = S.data.xsim;            
            S.data.mhempc.performance.ysim{num_sims} = S.data.ysim;
            S.data.mhempc.performance.xest{num_sims} = S.data.xest;
            S.data.mhempc.performance.ctrl{num_sims} = S.mpc.Controls;
            S.data.mhempc.performance.umeas{num_sims} = S.sensors.velocities;
            S.data.mhempc.performance.mpcRefs{num_sims} = S.data.references;
            S.data.mhempc.performance.exec_times{num_sims} = S.exec_time;
            S.data.mhempc.performance.obsPos{num_sims} = S.path.posObs;
            if S.config.save_workspace
                dirBase = '/home/nahuel/Documents/Nahuel/Algorithms/Obstacle-Avidance N-trailers/Simulations/sims Ts50ms/';
                nameFile = strcat([dirBase,'obsAvoidanceGNT-Nt',num2str(S.config.Nt),'-Nc',num2str(S.config.Nc),'-Np',num2str(S.config.Np),'-numMovObs',num2str(S.config.maxMovingObs),'-numStatObs',num2str(S.config.maxStaticObs),'-mthd-',S.config.method]);
                save(nameFile);
                %
                if S.config.lidarOn == true
                    stop(S.vlp16);
                    S = rmfield(S,'vlp16');
                    % S = rmfield(S,'obs');
                end
            end
        end
    end
end

function S = build_setup(S)
    % SIMULATION PARAMETERS ===============================================
    % S.ROS.IMU_ZEROING_VEC1  = -2.68;
    % S.ROS.IMU_ZEROING_VEC2  = -0.46;
    % S.ROS.IMU_ZEROING_MIC1  = 2.11;
    S.ROS.ENCODERS_PULSE_TO_RAD = 0.6 * pi / 180;
    S.config.ZEROING_VEC1   = true;
    S.config.ZEROING_VEC2   = false;
    S.config.ZEROING_MIC1   = false;
    S.config.num_meas_zeroing = 100;
    % Mobile Sensors ______________________________________________________
    S.mobile.useMobileSensors = true;
    S.mobile.nroMobileSensors = 1;
    % Solver ______________________________________________________________
    S.config.solver         = 'casadi'; % options: 'casadi', 'acado'    
    S.config.method         = 'proposed'; % 'proposed', 'fnmppc'
    S.config.integrator     = 'RK4'; % euler
%     if strcmp(S.config.method,'fnmppc')
%         S.config.integrator     = 'euler'; % euler
%     end
    % Simulation or field experiment (Husky)_______________________________
    S.config.SIM            = true;
    % Compute the Up-To-N-Steps Controllable set
    S.config.stab.compFlg   = false; % "Smplified design of practically stable MPC schemes", R. Comelli et. al.
    S.config.stab.nx0Stab   = 5000;
    % X_{N}^{\Omega_I}
    S.config.stab.betaMin   = -110*pi/180;
    S.config.stab.betaMax   = 110*pi/180;
    S.config.stab.thetaMin  = -pi/10;
    S.config.stab.thetaMax  = pi/10;
    S.config.stab.xMin      = -0.25;
    S.config.stab.xMax      = 0.25;
    S.config.stab.yMin      = -0.25;
    S.config.stab.yMax      = 0.25; 
    S.config.Omega_O        = Polyhedron('lb',[S.config.stab.xMin;S.config.stab.yMin],'ub',[S.config.stab.xMax;S.config.stab.yMax]);
    S.config.Omega_O_updt   = S.config.Omega_O;
    % X
    S.config.stab.XxMin     = -1;
    S.config.stab.XxMax     = 1;
    S.config.stab.XyMin     = -1;
    S.config.stab.XyMax     = 1;
    S.config.stab.XthetaMin = -pi;
    S.config.stab.XthetaMax = pi;
    % From numerical computations, I foun that for Nc=15, the currrent
    % control's costraints, a 2D projection of $X_N{\Omega_i}$ for beta_1
    % in [-1.2 1.2] (rad) is:
%     V                       = [-2 0; -1.6364 -0.1818; -1.6364 0.1818; 2 -0.1818; 2 0.1818];

    %
    S.config.stab.setNc     = 15;%[1 5 10 15 20 25 30 35 40];
    %
    S.config.Nt             = 1;
    V                       = [-5 -5; -5 5; 5 5; 5 -5]./(1+log10(S.config.Nt));
    S.config.X_N_Omega_I    = Polyhedron(V);
    S.config.X_N_Omega_I_updt = S.config.X_N_Omega_I;
    S.config.verbose        = false;
    S.config.calcMtxConvCoord = true; if ~S.config.calcMtxConvCoord; warning('WARNING!!! The matrix for correcting x-y coordinates is not being computed...'); end;
    S.config.t0             = 0;
    S.config.tf             = 3500;
    S.config.Ts             = 0.05;
    S.config.same_seed      = false;
    S.config.EXPORT         = false; if ~S.config.EXPORT; warning('WARNING!!! MHE and MPC were not compiled...'); end
    S.config.iters          = 0;
    S.config.time           = 0;
    S.config.NUM_SIMS       = length(S.config.stab.setNc);
    S.config.updtAC         = false;
    S.config.outputs        = [1:S.config.Nt+1,S.config.Nt+2:S.config.Nt+3];   % 2*S.config.Nt+6:2*S.config.Nt+7];
    S.config.Nc             = S.config.stab.setNc(S.data.num_sims);%6;                                               % Lenght of the control horizon
    S.config.Np             = 65;  % prediction horizon
    S.config.Ne             = 6;                                                % Lenght of the estimation window
    S.config.segmentTosteer = 0;

    S.config.maxMovingObs   = 0;
    S.config.maxStaticObs   = 0;
    S.config.totNumObs      = S.config.maxStaticObs+S.config.maxMovingObs;
    
    S.config.obsRand        = false;
    S.config.numStObsSgmnt  = 3;    % number of points to be evaluated in the parametric equation of the line
    S.config.numMvObsSgmnt  = 2;
    S.config.maxDistStatObs = 30;    % Static osbtacles beyond this distance are discarded
    S.config.numMeasEll     = 4;
    
    S.config.nroFrDetMObs   = S.config.numMeasEll;
    S.config.thrMovingObs   = 0.002;%0.008;
    S.config.distToUpdEll   = 0.2;
    S.config.accelOneMovObs = true;
    % Estimator and Control algorithm to execute ________________________
    S.config.obs_strategy   = 'gauss';
    S.config.mpc            = true;
    S.config.mhe            = true;
    % Disturbances ______________________________________________________
    S.config.noise          = 'gaussian'; % 'gaussian' or 'uniform'
    S.config.noise_lvl      = 0.*[(0.2*pi/180).*ones(S.config.Nt, 1); 0.2*pi/180; [0.05; 0.05]];%; [0.05; 0.05]]; % measurement noise amplitude: [betas'; theta_0; x_0; y_0; w0; v0]
    S.config.initUncertainty = 10*pi/180;
    S.config.slip           = [1; 1];
    S.config.procDist.type  = 'normal';
    S.config.procDist.amp   = 0.0.*[(0.5*pi/180).*ones(S.config.Nt,1); zeros(S.config.Nt+1,1); zeros(2,1); zeros(2,1)];
    S.config.gpsOutlier     = false;
    S.config.gpsOutlierAmp  = 2;
    S.config.gpsOutlierProbab = 10;
    S.config.IMUbias        = false;
    S.config.IMUbiasDev     = 5/180;
    S.config.model.uncty    = false;
    S.config.model.dev      = 55;                                       % porcentual value of uncertainty
    S.config.iterCorrecOn   = 200;
    S.config.CtrlNoise_lvl  = [0; 0];
    % Reference velocities ________________________________________________
    S.config.vNr            = 0;
    S.config.vNTarget       = 0.3;
    S.config.timeToReachVNr = 5; % seconds
    S.config.deltavN        = S.config.vNTarget/(S.config.timeToReachVNr/S.config.Ts + 2*(S.config.Ne+1));
    S.config.disLastTgt     = 0.75;
    % GPS's relative position to centre of vehicle ________________________
    S.config.gps_x          = -0.3;
    S.config.gps_y          = 0;
    S.config.gps_d          = sqrt(S.config.gps_x^2 + S.config.gps_y^2);
    S.config.gps_fcx        = 0;
    S.config.gps_fcy        = 0;
    S.config.SIMGPS         = false; % use to simulate date from GPS when it does not have signal. Just for indoor testing purposes.
    % Obstacles ___________________________________________________________
    S.config.lidarLimits.X  = 5;%1.0;       % distance to points from the lidar
    S.config.lidarLimits.Y  = 10;%2.3;   
    S.config.lidarLimits.Z  = 1;
    % boundaries for ignoring obstacles
    S.config.x_max          = 13;
    S.config.x_min          = 2;
    S.config.y_max          = 8;
    S.config.y_min          = 1.5;
    S.config.lidarOn        = false;
    S.config.obsDetection   = false;    
    % Plot data online during experiments. USeful for debbuging purposes __
    S.config.dataOnline     = true;
    S.config.save_workspace = false;    
    S.config.video          = false;
    %    
end

function S = setNc(S)
    S.config.Nc = S.config.stab.setNc(S.data.num_sims);
end

function S = gen_InitConditionsStab(S)    
    S.stability.xRange      = (S.config.stab.XxMin:0.1:S.config.stab.XxMax);
    S.stability.yRange      = (S.config.stab.XyMin:0.1:S.config.stab.XyMax);
    S.stability.theta0Range = (S.config.stab.XthetaMin:5*pi/180:S.config.stab.XthetaMax);
    S.stability.theta1Range = [-pi/2 -pi/4 0 pi/4 pi/2];
    S.stability.theta2Range = [-pi/2 -pi/4 0 pi/4 pi/2];
    indx        = 1;
    for x=1:length(S.stability.xRange)
        for y=1:length(S.stability.yRange)
            for theta0=1:length(S.stability.theta0Range)
                for theta1=1:length(S.stability.theta1Range)
%                     for theta2=1:length(S.stability.theta2Range)
%                         S = gen_x0(S,[S.stability.xRange(x),S.stability.yRange(y),S.stability.theta0Range(theta0),S.stability.theta0Range(theta0)+S.stability.theta1Range(theta1),S.stability.theta0Range(theta0)+S.stability.theta1Range(theta1)+S.stability.theta2Range(theta2)]);
                        S = gen_x0(S,[S.stability.xRange(x),S.stability.yRange(y),S.stability.theta0Range(theta0),S.stability.theta0Range(theta0)+S.stability.theta1Range(theta1)]);
                        S.initConditionsStab{indx} = S.init_condition.x0;
                        indx = indx+1;
%                     end
                end
            end
        end        
    end  
    %
%     S.initConditionsStab = gpuArray(S.initConditionsStab);
    S.config.stab.nx0Stab = indx-1;
end

function S = solveMPCStab(S)
%     setReference(S.mpc.mpcCasadi,[0;0;0]);
tic
mpc = S.mpc.mpcCasadi;
X0 = S.initConditionsStab;
    parfor i=1:S.config.stab.nx0Stab
%         S = init_mpc(S);
        %
        setReference(mpc,[0;0;0]);
        setU_l(mpc,[0;0]);
        setq0(mpc,X0{i});
        solve(mpc);
%         S.trajStab{(S.data.num_sims-1)*S.config.stab.nx0Stab + i}   = S.mpc.mpcCasadi.Qtraj;
%         S.UStab{(S.data.num_sims-1)*S.config.stab.nx0Stab + i}   = S.mpc.mpcCasadi.Utraj;
%         S.xNcStab       = [S.xNcStab, S.mpc.mpcCasadi.Qtraj(:,end)];
        %
        S.config.stab.nx0Stab-i
    end
end

function S = CALL_SIMULATOR(S)
    % Solve estimation problem ________________________________
    S = call_ESTIMATOR(S);
    % Obstacle detection ______________________________________
    S = call_OBSDETECTOR(S);
    % Path-tracking problem ___________________________________
    S = call_PFA(S);
    % Solve control problem ___________________________________
    S = call_CONTROLLER(S);                
    % Apply controls to the Husky _____________________________
    S = call_SYSTEM(S);                
    % Update measurements _____________________________________
    S = update_measurements(S);
    % Perform real-time iteration _____________________________
    S = call_RTF(S);
end

function S = call_RTF(S)
    tic;
    S.exec_time.t_tot   = [ S.exec_time.t_tot, S.exec_time.t_mhe(end)+S.exec_time.t_pfa(end)+S.exec_time.t_mpc(end)+...
                            S.exec_time.t_ctrl(end)+S.exec_time.t_mis(end)+S.exec_time.t_obsdetector(end)+...
                            S.exec_time.t_sensors(end)+S.exec_time.t_plt(end)];
    S.exec_time.t_acum  = S.exec_time.t_tot(end);
    while ((S.exec_time.t_acum + toc) < S.config.Ts) && ~S.config.SIM; end
end

function plot_mono2(S, qk, clr)
    if nargin == 3
        clrTractor          = clr;
        clrWheel            = clr;
        clrAxe              = clr;
        clrLongAxe          = clr;
        clrTrailerLoad      = clr;
        clrLastTrailerLoad  = clr;
        noRef               = true;
    else
        clrTractor          = 'y';
        clrWheel            = 'k';
        clrAxe              = 'k';
        clrLongAxe          = 'k';
        clrTrailerLoad      = 'b';
        clrLastTrailerLoad  = 'r';
        noRef               = false;
    end
    % Declare all local variables for paralellising
    if strcmp(S.config.method,'proposed')
        if ~isempty(S.mpc.mpcCasadi.Qtraj)
            Qtraj = S.mpc.mpcCasadi.Qtraj;
        else
            Qtraj = repmat(S.mpc.mpcCasadi.q0bar, 1, S.config.Np);
        end
    end
    if strcmp(S.config.method,'fnmppc')
        if ~isempty(S.fnmppc.Qtraj)
            Qtraj = S.fnmppc.Qtraj;
        else
            Qtraj = repmat(S.fnmppc.optiMPC.value(S.fnmppc.qinit), 1, S.config.Np);
        end
    end    
    coordinates         = S.path.coordinates;
    Nt                  = S.config.Nt;
    listOfObs           = S.path.listOfObs;
    safeMargin          = S.path.safeMargin;
    vehicleDims         = S.path.vehicleDims;
    dynamicObs          = S.path.dynamicObs;
    segmentTosteer      = S.config.segmentTosteer;
    ref                 = S.controller.ref;
    xReachable          = ref.xReachable;
    yReachable          = ref.yReachable;
    nearObs             = S.path.nearObs;
    THETA               = S.mpc.mpcCasadi.optiMPC.value(S.mpc.mpcCasadi.THETA);
    XC                  = S.mpc.mpcCasadi.optiMPC.value(S.mpc.mpcCasadi.XC);
    YC                  = S.mpc.mpcCasadi.optiMPC.value(S.mpc.mpcCasadi.YC);
    A                   = S.mpc.mpcCasadi.optiMPC.value(S.mpc.mpcCasadi.A);
    B                   = S.mpc.mpcCasadi.optiMPC.value(S.mpc.mpcCasadi.B);
    T0                  = S.mpc.mpcCasadi.optiMPC.value(S.mpc.mpcCasadi.T0);
    KT                  = S.mpc.mpcCasadi.optiMPC.value(S.mpc.mpcCasadi.KT);
    flg                 = S.path.occluded.flg;
    p1                  = S.path.occluded.p1;
    p2                  = S.path.occluded.p2;
    p3                  = S.path.occluded.p3;
    XYtrailerAxe        = S.system.XYtrailerAxe;
    XYtrailerWheelLeft  = S.system.XYtrailerWheelLeft;
    XYtrailerWheelRight = S.system.XYtrailerWheelRight;
    XYtrailerLongAxe    = S.system.XYtrailerLongAxe;
    

    if ~noRef        
        % Plot the path
        plot(coordinates(1,:),coordinates(2,:),'color',[0.7 0.7 0.7],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
        % Plot the predicted trajectory for each segment of the N-trailer
        for i=1:Nt+1
            if i==1
                plot(Qtraj(2*Nt+1+(i-1)*2+1,:),Qtraj(2*Nt+1+i*2,:),'y','LineWidth',4);
            elseif i~=Nt+1
                plot(Qtraj(2*Nt+1+(i-1)*2+1,:),Qtraj(2*Nt+1+i*2,:),'b','LineWidth',4);
            else
                plot(Qtraj(2*Nt+1+(i-1)*2+1,:),Qtraj(2*Nt+1+i*2,:),'r','LineWidth',4);
            end
        end
        % Plot the static obstacles as outer circles
        if~isempty(listOfObs)
            nroObs = size(listOfObs,1);
            for i=1:nroObs
                outter1 = circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3),'color',[6 18 131]./255);
                alpha(outter1,0.2)
                outter2 = circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3)-safeMargin,'color',[19 141 144]./255);
                alpha(outter2,0.4)
                circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3)-safeMargin-vehicleDims,'color',[253 60 60]./255);
            end
        end
        % Plot the moving obstacles as outer circles
        if~isempty(dynamicObs)
            nroObs = size(dynamicObs,1);
            for i=1:nroObs
                outter1 = circles(dynamicObs(i,1),dynamicObs(i,2),dynamicObs(i,3),'color',[0 207 250]./255);
                alpha(outter1,0.2);
                outter2 = circles(dynamicObs(i,1),dynamicObs(i,2),dynamicObs(i,3)-safeMargin,'color',[76 63 84]./255);
                alpha(outter2,0.4);
                circles(dynamicObs(i,1),dynamicObs(i,2),dynamicObs(i,3)-safeMargin-vehicleDims,'color',[255 0 56]./255);
            end
        end
    end    
    % Plot the current reference
    if segmentTosteer == 0
        clr         = 'y';
        clrReach    = [0.95 0.95 0.1];
    elseif segmentTosteer == Nt
        clr         = 'r';
        clrReach    = [1 0.1 0.1];
    else
        clr = 'b';
        clrReach    = [0.1 0.1 1];
    end
    plot(ref.x,ref.y,clr,'marker','+','linewidth',2,'markersize',15)
    plot(ref.x,ref.y,clr,'marker','o','linewidth',2,'markersize',15)
    % plot the reachable reference
    plot(xReachable,yReachable,'color',clrReach,'marker','+','linewidth',2,'markersize',15)
    plot(xReachable,yReachable,'color',clrReach,'marker','o','linewidth',2,'markersize',15)
    %
    betas   = qk(1:Nt);    
    thetas  = qk(Nt+1:2*Nt+1); 
    xy0     = qk(2*Nt+2:2*Nt+3);
    xyi     = zeros(2,Nt+1);
    xyi(:,1) = xy0;
    for i=1:Nt
        xy_aux      = qk(2*Nt+3+(i-1)*2+1:2*Nt+3+(i)*2);
        xyi(:,i+1)  = xy_aux;
        %
        pltLine = [];
        if ~isempty(nearObs) % Line betwwen trauiler and near obstacles
%             for j=1:size(nearObs{i+1},1)
%                 xy      = nearObs{i+1}(j,1:2)';
%                 pltLine = [pltLine, [xy_aux, xy]];
%             end
            if ~isempty(pltLine)
%                 if i==Nt
%                     line(pltLine(1,:),pltLine(2,:),'color',[1 0 0]);
%                     alpha(0.4)
%                 else
%                     line(pltLine(1,:),pltLine(2,:),'color',[0 0 1]);
%                 end
            end
        end        
    end
%     if ~isempty(nearObs)%{1}
%         pltLine = [];
%         for j=1:size(nearObs{1},1)
%             xy      = nearObs{1}(j,1:2)';
%             pltLine = [pltLine, [xy0, xy]];
%         end
%         if ~isempty(pltLine)
%             line(pltLine(1,:),pltLine(2,:),'color',[1 1 0]);
%         end
%     end
    % Plot moving obstacle
    t1      = 0:0.1:2*pi;    
    for i=1:length(A)
%         x = cos(THETA(i))*XC(i)+A(i)*cos(THETA(i)).*cos(t1)-sin(THETA(i))*YC(i)-B(i)*sin(THETA(i)).*sin(t1);
%         y = sin(THETA(i))*XC(i)+A(i)*sin(THETA(i)).*cos(t1)+cos(THETA(i))*YC(i)+B(i)*cos(THETA(i)).*sin(t1);
%         plot(x,y,'color',[102 165 173]./255);
        if ~isempty(KT)
            t = 0:S.config.Np-1;
            x = XC(i) + A(i)*cos(T0(i)+KT(i).*t);
            y = YC(i) + B(i)*sin(T0(i)+KT(i).*t);
            plot(x,y,'color',[161 214 226]./255,'LineWidth',2);
        end
    end
    %
    % plot occlusion's triangle
    if flg
        T = delaunayTriangulation([p1;p2;p3]);
        if ~isempty(T.ConnectivityList)
            hold on;
            triplot(T,'color',[72 151 216]./255);
            hold off;
        end
    end 
    % Plot set $\Omega_I=\Omega_O$
%     hold on;
% %     R           = [cos(-ref.theta) -sin(-ref.theta); sin(-ref.theta) cos(-ref.theta)];
% %     Omega_O     = S.config.Omega_O_updt;%*R + [ref.x; ref.y];
%     O_O         = S.config.Omega_O_updt.plot('color','r'); 
%     alpha(O_O,0.1);    
%     % Plot set $X_N^{\Omega_I}$
% %     R           = [cos(-thetas(1)) -sin(-thetas(1)); sin(-thetas(1)) cos(-thetas(1))];
% %     X_N_Omega_I = S.config.X_N_Omega_I_updt;%*R + xy0;
%     X_N         = S.config.X_N_Omega_I_updt.plot('color','b');
%     alpha(X_N,0.05);
%     hold off;
    %
    R                       = [cos(thetas(1)-pi/2), -sin(thetas(1)-pi/2); sin(thetas(1)-pi/2), cos(thetas(1)-pi/2)];
    Rh                      = [cos(-thetas(1)-pi/2), -sin(-thetas(1)-pi/2); sin(-thetas(1)-pi/2), cos(-thetas(1)-pi/2)];

    tractorBodyPlt          = S.system.XYtracBody*Rh + xyi(:,1);
    hold on;
    tractorBodyPlt.plot('color','y');
    hold off;


    tractorWheelLeftPlt     = R * S.system.XYtracWheelLeft + xyi(:,1);
    tractorWheelRightPlt    = R * S.system.XYtracWheelRight + xyi(:,1);
    tractorAxePlt           = R * S.system.XYtracAxe + xyi(:,1);
    %
%     line(tractorBodyPlt(1,:),tractorBodyPlt(2,:),'color',clrTractor,'linewidth',3);
    line(tractorWheelLeftPlt(1,:),tractorWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorWheelRightPlt(1,:),tractorWheelRightPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorAxePlt(1,:),tractorAxePlt(2,:),'color',clrAxe,'linewidth',2);
    %
    for i=1:Nt
        R                       = [cos(thetas(i+1)-pi/2), -sin(thetas(i+1)-pi/2); sin(thetas(i+1)-pi/2), cos(thetas(i+1)-pi/2)];
        Rtr                     = [cos(-thetas(i+1)+pi/2), -sin(-thetas(i+1)+pi/2); sin(-thetas(i+1)+pi/2), cos(-thetas(i+1)+pi/2)];
        trailerAxePlt           = R * XYtrailerAxe + xyi(:,i+1);
        trailerWheelLeftPlt     = R * XYtrailerWheelLeft + xyi(:,i+1);
        trailerWheelRightPlt    = R * XYtrailerWheelRight + xyi(:,i+1);
        trailerLongAxePlt       = R * XYtrailerLongAxe((i-1)*2+1:i*2,:) + xyi(:,i+1);
%         trailerLoadPlt          = R * S.system.XYtrailerLoad((i-1)*2+1:i*2,:) + xyi(:,i+1);
        trailerLoadPlt          = S.system.XYtrailerLoad{i}*Rtr + xyi(:,i+1);
        hold on;
        if i==Nt
            trailerLoadPlt.plot('color','r');
        else
            trailerLoadPlt.plot('color','b');
        end

        if i~= Nt
            trailerClr = clrTrailerLoad;
        else
            trailerClr = clrLastTrailerLoad;
        end
        circles(trailerLongAxePlt(1,2),trailerLongAxePlt(2,2),S.system.r/2,'edgecolor',clrLongAxe,'facecolor','none')
        line([trailerLongAxePlt(1,2),xyi(1,i)],[trailerLongAxePlt(2,2),xyi(2,i)],'color',clrLongAxe,'linewidth',1)    
        line(trailerAxePlt(1,:),trailerAxePlt(2,:),'color',clrAxe,'linewidth',2)
        line(trailerWheelLeftPlt(1,:),trailerWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerWheelRightPlt(1,:),trailerWheelRightPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerLongAxePlt(1,:),trailerLongAxePlt(2,:),'color',clrLongAxe,'linewidth',2)
%         line(trailerLoadPlt(1,:),trailerLoadPlt(2,:),'color',trailerClr,'linewidth',2)
    end
    %
    if ~noRef
        hold off;
    end    
    %
    xlim([min(coordinates(1,:))-2 max(coordinates(1,:))+2]); ylim([min(coordinates(2,:))-4 max(coordinates(2,:))+4]);
    %
    drawnow limitrate    
end

function S = call_init_functions(S)
    S = reserve_temp_memory(S);
    S = gen_init_conditions(S);
    S = init_flags_and_counters(S);    

    % temporarily here:
    if S.config.lidarOn
        S = init_obsDetection(S);
        start(S.vlp16);
    end
end

function S = compute_tras_rot_mtx(S)
    sdpvar a11 a12 a21 a22;
    sdpvar x1bar y1bar x2bar y2bar;
    % Besides the reference point, two more are needed to obtaint the local
    % reference frame
    % AZOTEA AC3E *********************************************************
    % Coordinates of point (x, 0)
    S.ROS.local_coord_1.lat     = -33.034213;       % hand coded value, measurement from rtk
    S.ROS.local_coord_1.long    = -71.592168;   % hand coded value, measurement from rtk
    % CANCHA DE FUTBOL ****************************************************
%     S.ROS.local_coord_1.lat     = -33.03517;       % hand coded value, measurement from rtk
%     S.ROS.local_coord_1.long    = -71.594195;   % hand coded value, measurement from rtk
    % *********************************************************************
    [x1, y1]                    = latlon2xy(S.ROS.local_coord_1.lat, S.ROS.local_coord_1.long, S.ROS.LAT0, S.ROS.LON0);
    x1                          = x1*1000;
    y1                          = y1*1000;
    x                           = norm([x1 y1]);
    % AZOTEA AC3E *********************************************************
    % Coordinates of point (0, y)
    S.ROS.local_coord_2.lat     = -33.034088;        % hand coded value, measurement from rtk
    S.ROS.local_coord_2.long    = -71.59211333;         % hand coded value, measurement from rtk
    % CANCHA DE FUTBOL ****************************************************
%     S.ROS.local_coord_2.lat     = -33.034403333;        % hand coded value, measurement from rtk
%     S.ROS.local_coord_2.long    = -71.594075;         % hand coded value, measurement from rtk
    % *********************************************************************
    [x2, y2]                    = latlon2xy(S.ROS.local_coord_2.lat, S.ROS.local_coord_2.long, S.ROS.LAT0, S.ROS.LON0);
    x2                          = x2*1000;
    y2                          = y2*1000;
    y                           = norm([x2 y2]);
    % With the "origin" and a point alongside each axe, compute the mtx
    A                           = [a11 a12; a21 a22];
    v                           = [x1; y1; x2; y2];
    b                           = [x1bar; y1bar; x2bar; y2bar];
    Constraints                 = [[A zeros(2); zeros(2) A]*v - b == zeros(4,1); x1bar*x2bar + y1bar*y2bar == 0];
    %
    obj                         = (x1bar - x)^2 + (y2bar - y)^2;
    % 
    optimize(Constraints, obj);
    %
    S.ROS.Mtx                   = value(A);
end

function S = update_measurements(S)
    tic;
    if S.config.SIM == true
        if strcmp(S.config.noise,'gaussian')
            noise = randn(S.system.ny,1).*S.config.noise_lvl;
        else
            noise = (2*rand(S.system.ny,1)-1).*S.config.noise_lvl;
        end
        if S.config.gpsOutlier == true
            if randi(S.config.gpsOutlierProbab) == 1
                noise(S.config.Nt+2:S.config.Nt+3) = S.config.gpsOutlierAmp*randn(2,1);
            end
        end
        if S.config.IMUbias == true
            noise(S.config.Nt+1) = noise(S.config.Nt+1) + S.config.IMUbiasDev*S.config.iters*S.config.Ts/3600; % deviation  per hour     %*S.data.xsim(S.config.Nt+1,end);
        end

        S.data.measNoise = [S.data.measNoise, noise];
        S.data.ysim = [S.data.ysim, S.data.xsim(S.config.outputs,end) + noise];
        %
        unoise = randn(S.system.nu,1).*S.config.CtrlNoise_lvl;
        S.data.UmeasNoise = [S.data.UmeasNoise, unoise];
        S.sensors.velocities = S.mpc.Controls(:,end) + unoise;
    else
        S           = read_sensors(S);
        S           = measurements_vector(S);
        S.data.ysim = [S.data.ysim, S.sensors.measurement]; % <-- MEASUREMENTS are stored in this vector in this order: [beta_1, beta_2, theta_0, x_0, y_0]
    end
    S = update_EstimatorMeasurement(S);
    S.exec_time.t_sensors = [S.exec_time.t_sensors, toc];
end

function flag = check_end_condition(S)
    flag        = false;
    condition = norm(S.data.xest(2*S.config.Nt+2:2*S.config.Nt+3,end)-S.path.coordinates(:,end));
    if (condition <= S.config.disLastTgt) && S.path.last_tgt
        flag = true;
    end
end

function S = ROS(S)
    S = init_ROS(S);
    S = create_obj_sens(S);
    S = create_obj_vel(S);    
    S = zeroing_imu(S);
end

function S = MOBILE(S)
    S = init_mobile(S);
end

function S = init_mobile(S)
    conn_devices = size(mobiledevlist,1);
    while conn_devices < S.mobile.nroMobileSensors
        conn_devices = size(mobiledevlist,1);
        fprintf('Waiting for mobiles to connect...\n')
    end
    S.mobile.devs    = {};
    S.mobile.devs{1} = 'POCO M4 Pro 5G';
    S.mobile.devs{2} = 'POCO M4 Pro 5G';
    S.mobile.devs{3} = 'POCO M4 Pro 5G';
    S.mobile.devs{4} = 'POCO M4 Pro 5G';
    %
    S.mobile.sensors = {};
    for i=1:S.config.nroMobileSensors
        S.mobile.sensors{i} = mobiledev(S.mobile.devs{i});
    end
      
end

function S = readMobile(S)
    for i=1:S.config.nroMobileSensors
        data = [S.mobile.sensors{i}.Acceleration'; S.mobile.sensors{i}.AngularVelocity';...
            S.mobile.sensors{i}.MagneticField'; S.mobile.sensors{i}.Orientation';...
            S.mobile.sensors{i}.Latitude; S.mobile.sensors{i}.Longitude; S.mobile.sensors{i}.Speed;...
            S.mobile.sensors{i}.Course; S.mobile.sensors{i}.Altitude; S.mobile.sensors{i}.HorizontalAccuracy ];
        S.mobile.data{i} = [S.mobile.data{i}, data];
    end
    
end

function S = init_ROS(S)
    % Init coordinates of my local 2D plane
    S = get_reference(S);
    
    if S.config.calcMtxConvCoord == true
        S = compute_tras_rot_mtx(S);
    else
        S.ROS.Mtx = eye(2);
    end
    
    S.sensors.measurements = [];
    % Create and init Husky's velocities
    S.mpc.Husky_v0 = 0;
    S.mpc.Husky_w0 = 0;
    % Here some instructions need to be followed from matlab's terminal
    rosshutdown;
    input('use el comando roscore desde terminal...\n')
    %
    rosshutdown;
    rosinit;
    %
    fprintf('\n---------------------------------------------------------\n');
    fprintf('RECORDAR RESETEAR LOS ENCODERS ANTES DEL EXPERIMENTOS...\n');
    fprintf('---------------------------------------------------------\n');
    %
    fprintf(['\n\nroslaunch husky_base base.launch...\n' ...
        'roslaunch parches rtk.launch...\n' ...
        'roslaunch vectornav vectornav.launch...\n' ...
        'roslaunch parches parche_imu_ins.launch...\n' ...
        'roslaunch vectornav vectornav2.launch...\n' ...
        'roslaunch imu_3dm_gx4 imu.launch...\n' ...
        'roslaunch parches encoders.launch...\n' ...
        'roslaunch parches parche_speed_holder.launch...\n'])
    aux1 = input('Presione enter...\n');
end

function S = create_obj_sens(S)
    sub_rtk         = rossubscriber('/fix');
    sub_vec1        = rossubscriber('/vectornav/IMU');
    sub_vec1_vel    = rossubscriber('/imu_INS');
    % sub_vec2        = rossubscriber('/vectornav2/IMU');
    sub_encoders    = rossubscriber('/enc');
    % sub_micro1      = rossubscriber('/imu/pose');
%     sub_lidar       = rossubscriber('/angles');
    %
    S.ROS.rtk            = sub_rtk;
    S.ROS.vectornav1     = sub_vec1;
    % S.ROS.vectornav2     = sub_vec2;
    S.ROS.vectornav1_vel = sub_vec1_vel;
    % S.ROS.microstrain    = sub_micro1;
    S.ROS.IncEncoders    = sub_encoders;
%     S.ROS.betasFromLidar = sub_lidar;

    % COrection factor for unwraping phase on real-time
    S.ROS.pose.correction_vec1              = 0;      % This value should be substracted to every pose measurement
    S.ROS.sensors.vectornav_euler_vec1      = [];
    S.ROS.sensors.vectornav_euler_vec1_Old  = 0;     % Value for computinng the difference. If it is bigger in absolute value than pi, a correction is needed
    %
    % S.ROS.pose.correction_vec2              = 0;
    % S.ROS.sensors.vectornav_euler_vec2      = [];
    % S.ROS.sensors.vectornav_euler_vec2_Old  = 0;
    %
    % S.ROS.pose.correction_micro1            = 0;      % This value should be substracted to every pose measurement
    % S.ROS.sensors.microstrain_euler_micro1    = [];
    % S.ROS.sensors.microstrain_euler_micro1_Old  = 0;
end

function S = create_obj_vel(S)
    [pub,msg]       = rospublisher('/cmd_vel_aux','geometry_msgs/Twist');
    %
    S.ROS.CTRL.pub       = pub;
    S.ROS.CTRL.msg       = msg;
end

function S = get_reference(S)
    % Coordiantes of the origin of our local reference 2D-axes x-y. Values
    % charged by hand before carry out the experiments.
    % AZOTEA AC3E *********************************************************
    S.ROS.LAT0   = -33.034115;
    S.ROS.LON0  = -71.592205;
    % CANCHA DE FUTBOL ****************************************************
%     S.ROS.LAT0   = -33.03453;
%     S.ROS.LON0  = -71.594498;
    %
end

function S = read_rtk(S)
    if ~S.config.SIMGPS
        S.ROS.sensors.gpsrtk = receive(S.ROS.rtk);
        
        if isnan(S.ROS.sensors.gpsrtk.Latitude) || isnan(S.ROS.sensors.gpsrtk.Longitude)
            S = write_ctrl(S, 0, 0);
            while isnan(S.ROS.sensors.gpsrtk.Latitude) || isnan(S.ROS.sensors.gpsrtk.Longitude)
                fprintf('NO GPS SIGNAL...\n')
                pause(0.2);
            end
        end
        % Convert lat &v lon to x-y coordinated in a tangential plane ro
        % earth's surface (tangent to my reference point)
        [x,y]   = latlon2xy(S.ROS.sensors.gpsrtk.Latitude,S.ROS.sensors.gpsrtk.Longitude,S.ROS.LAT0,S.ROS.LON0);
        % Convert Km to m
        x_raw   = x*1000;
        y_raw   = y*1000;
        % Adjust to my local 2D plane through the rotation matrix
        xy_cor  = S.ROS.Mtx * [x_raw; y_raw];
        % The, store values to be used later
        S.ROS.sensors.rtk.x0 = xy_cor(1) - S.config.gps_fcx;
        S.ROS.sensors.rtk.y0 = xy_cor(2) - S.config.gps_fcy;
% pos=[S.ROS.sensors.rtk.x0, S.ROS.sensors.rtk.y0]
    else
        S.ROS.sensors.rtk.x0 = S.data.xsim(2*S.config.Nt+2,end);
        S.ROS.sensors.rtk.y0 = S.data.xsim(2*S.config.Nt+3,end);
        %
        fprintf('Simulated data from GPS...\n')
    end
end

function S = write_ctrl(S, w0, v0)
    tic;
    S.ROS.CTRL.msg.Linear.X    = v0;
    S.ROS.CTRL.msg.Angular.Z   = w0;
    send(S.ROS.CTRL.pub,S.ROS.CTRL.msg);
    S.exec_time.t_ctrl = [S.exec_time.t_ctrl, toc];
end

function S = zeroing_imu(S)
%
    fprintf('\nZeroing imus...\n');
%
    if S.config.ZEROING_VEC1
        zero_vec1  = 0;
        for i=1:S.config.num_meas_zeroing
            vec1_orientation    = receive(S.ROS.vectornav1);
            quat_vec1           = vec1_orientation.Orientation;
            raw_vec1            = quat2eul([quat_vec1.X,quat_vec1.Y,quat_vec1.Z,quat_vec1.W]);
            zero_vec1           = zero_vec1 + raw_vec1(3);
        end
        S.ROS.IMU_ZEROING_VEC1 = zero_vec1 / S.config.num_meas_zeroing;
    else
        S.ROS.IMU_ZEROING_VEC1 = 2.48;
    end
    %
    if S.config.ZEROING_VEC2
        zero_vec2     = 0;
        for i=1:S.config.num_meas_zeroing
            vec2_orientation    = receive(S.ROS.vectornav2);
            quat_vec2           = vec2_orientation.Orientation;
            raw_vec2            = quat2eul([quat_vec2.X,quat_vec2.Y,quat_vec2.Z,quat_vec2.W]);
            zero_vec2           = zero_vec2 + raw_vec2(3);
        end
        S.ROS.IMU_ZEROING_VEC2 = zero_vec2 / S.config.num_meas_zeroing;
    else
        S.ROS.IMU_ZEROING_VEC2 = -2.93;
    end
    %
    if S.config.ZEROING_MIC1
        zero_mic1     = 0;
        for i=1:S.config.num_meas_zeroing
            micro1_orientation  = receive(S.ROS.microstrain);
            quat_mic1           = micro1_orientation.Pose.Orientation;
            raw_mic1            = quat2eul([quat_mic1.X,quat_mic1.Y,quat_mic1.Z,quat_mic1.W]);
            zero_mic1           = zero_mic1 + raw_mic1(3);
        end
        S.ROS.IMU_ZEROING_MIC1 = zero_mic1 / S.config.num_meas_zeroing;
    else
        S.ROS.IMU_ZEROING_MIC1 = -2.96;
    end
    %
    fprintf('Zeroing imus ready!\n\n');
end

function S = read_vectornav(S) % measure theta0 and the speeds
    % pose
    S.ROS.sensors.vectornav             = receive(S.ROS.vectornav1);
    quat                                = S.ROS.sensors.vectornav.Orientation;
    S.ROS.sensors.vectornav_euler_vec1  = quat2eul([quat.X,quat.Y,quat.Z,quat.W]);    
    % Unwrap phase from online data _______________________________________
    if (S.ROS.sensors.vectornav_euler_vec1(3)-S.ROS.sensors.vectornav_euler_vec1_Old) >= pi
        S.ROS.pose.correction_vec1 = S.ROS.pose.correction_vec1 + 2*pi;
    elseif (S.ROS.sensors.vectornav_euler_vec1(3)-S.ROS.sensors.vectornav_euler_vec1_Old) <= -pi
        S.ROS.pose.correction_vec1 = S.ROS.pose.correction_vec1 - 2*pi;         
    end         
    %
    S.ROS.sensors.vectornav_euler_vec1_Old = S.ROS.sensors.vectornav_euler_vec1(3);    
    % Do Not compute the attitude angle in my reference frame -------------
    S.ROS.sensors.vectornav_theta0 = -S.ROS.sensors.vectornav_euler_vec1(3) + S.ROS.IMU_ZEROING_VEC1 + S.ROS.pose.correction_vec1;
    % Measure the speed ---------------------------------------------------
    S.ROS.sensors.vectornav_vel     = receive(S.ROS.vectornav1_vel);
    S.ROS.sensors.vectornav_NedVelX = S.ROS.sensors.vectornav_vel.Data(1);
    S.ROS.sensors.vectornav_NedVelY = S.ROS.sensors.vectornav_vel.Data(2);
    S.ROS.sensors.vectornav_v0      = sqrt(S.ROS.sensors.vectornav_NedVelX^2 + S.ROS.sensors.vectornav_NedVelY^2);
    S.ROS.sensors.vectornav_w0      = -1 * S.ROS.sensors.vectornav.AngularVelocity.Z;
    % Compute correction factor ___________________________________________
    S.config.gps_fcx = S.config.gps_d * sin(S.ROS.sensors.vectornav_theta0);
    S.config.gps_fcy = S.config.gps_d * cos(S.ROS.sensors.vectornav_theta0);
end

function S = read_vectornav2(S) % measure theta2
    % pose
    S.ROS.sensors.vectornav2            = receive(S.ROS.vectornav2);
    quat                                = S.ROS.sensors.vectornav2.Orientation;
    S.ROS.sensors.vectornav_euler_vec2  = quat2eul([quat.X,quat.Y,quat.Z,quat.W]);    
    % Unwrap pahse from online data _______________________________________
    if (S.ROS.sensors.vectornav_euler_vec2(3)-S.ROS.sensors.vectornav_euler_vec2_Old) >= pi
        S.ROS.pose.correction_vec2 = S.ROS.pose.correction_vec2 + 2*pi;
    elseif (S.ROS.sensors.vectornav_euler_vec2(3)-S.ROS.sensors.vectornav_euler_vec2_Old) <= -pi
        S.ROS.pose.correction_vec2 = S.ROS.pose.correction_vec2 - 2*pi;         
    end         
    %
    S.ROS.sensors.vectornav_euler_vec2_Old = S.ROS.sensors.vectornav_euler_vec2(3);
    % No compute the attitude angle in my reference frame -----------------
    S.ROS.sensors.vectornav_theta2 = -S.ROS.sensors.vectornav_euler_vec2(3) + S.ROS.IMU_ZEROING_VEC2 + S.ROS.pose.correction_vec2;
end

function S = read_microstrain(S)
    S.ROS.sensors.microstrain   = receive(S.ROS.microstrain);
    quat                        = S.ROS.sensors.microstrain.Pose.Orientation;
    S.ROS.sensors.microstrain_euler_micro1 = quat2eul([quat.X,quat.Y,quat.Z,quat.W]);
    % Unwrap pahse from online data _______________________________________
    if (S.ROS.sensors.microstrain_euler_micro1(3)-S.ROS.sensors.microstrain_euler_micro1_Old) >= pi
        S.ROS.pose.correction_micro1 = S.ROS.pose.correction_micro1 + 2*pi;
    elseif (S.ROS.sensors.microstrain_euler_micro1(3)-S.ROS.sensors.microstrain_euler_micro1_Old) <= -pi
        S.ROS.pose.correction_micro1 = S.ROS.pose.correction_micro1 - 2*pi;         
    end         
    %
    S.ROS.sensors.microstrain_euler_micro1_Old = S.ROS.sensors.microstrain_euler_micro1(3);
    % No compute the attitude angle in my reference frame -----------------
    S.ROS.sensors.microstrain_theta1 = -S.ROS.sensors.microstrain_euler_micro1(3) + S.ROS.IMU_ZEROING_MIC1 + S.ROS.pose.correction_micro1;
end

function S = read_encoders(S)
    S.ROS.sensors.encoder.data      = receive(S.ROS.IncEncoders);    
    S.ROS.sensors.encoders.beta1    = double(S.ROS.sensors.encoder.data.Data(1)) * S.ROS.ENCODERS_PULSE_TO_RAD;
    S.ROS.sensors.encoders.beta2    = double(S.ROS.sensors.encoder.data.Data(2)) * S.ROS.ENCODERS_PULSE_TO_RAD;
end

function S = read_sensors(S)
    S = read_vectornav(S);
    % S = read_vectornav2(S);
    % S = read_microstrain(S);
    S = read_encoders(S);
    S = read_rtk(S);
end

function S = measurements_vector(S)
%     S.sensors.measurement   = [S.ROS.sensors.encoders.beta1; S.ROS.sensors.encoders.beta2; S.ROS.sensors.vectornav_theta0; S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0; S.ROS.sensors.vectornav_w0; S.ROS.sensors.vectornav_v0];
    S.sensors.measurement   = [S.ROS.sensors.encoders.beta1; S.ROS.sensors.encoders.beta2; S.ROS.sensors.vectornav_theta0; S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0];
    S.sensors.velocities    = [S.sensors.velocities, [S.ROS.sensors.vectornav_w0; S.ROS.sensors.vectornav_v0]];
    %
    S.sensors.measurements  = [S.sensors.measurements, S.sensors.measurement];   
    % Store the attitude of each trailer
    S.sensors.theta0        = [S.sensors.theta0, S.ROS.sensors.vectornav_theta0];
    % S.sensors.theta1        = [S.sensors.theta1, S.ROS.sensors.microstrain_theta1];
    % S.sensors.theta2        = [S.sensors.theta2, S.ROS.sensors.vectornav_theta2];
end

function S = update_EstimatorMeasurement(S)
%     if S.config.mhe
%         updateMeasurement(S.algorithms.mheCasadi, double(S.data.ysim(:,end)));    % in S.data.ysim(:,end-2) I save the velocities too
%         updateInput(S.algorithms.mheCasadi, double(S.sensors.velocities(:,end)));
%     else
%         updateMeasurement(S.algorithms.ekfCasadi, double(S.data.ysim(1:end-2,end)));
%         updateInput(S.algorithms.ekfCasadi, double(S.sensors.velocities(:,end)));
%     end
end

function S = call_ESTIMATOR(S)
    tic;
%     if S.config.mhe
%         solve(S.algorithms.mheCasadi);
%         q_k = S.algorithms.mheCasadi.q_k;
%     else
%         solve(S.algorithms.ekfCasadi);
%         q_k = S.algorithms.ekfCasadi.x_k;
%     end
q_k = S.data.xsim(:,end);
    S.data.xest = [S.data.xest, q_k];
    S.exec_time.t_mhe  = [S.exec_time.t_mhe, toc];
end

function S = updateTarget(S)
    disAcum = 0;
    kMax    = length(S.path.coordinates)-1;
    vDir    = S.path.coordinates(:,end)-S.path.coordinates(:,end-1);
    for i=S.path.acumIndx+1:kMax
        vDir = S.path.coordinates(:,i)-S.path.coordinates(:,i-1);
        disAcum = disAcum + norm(vDir);
        if (disAcum/S.config.Ts >= S.config.vNTarget*S.controller.modVelFac) || (i==kMax)
            S.path.acumIndx     = i;
            if i==kMax
                S.path.last_tgt = true;
            end
            break;
        end        
    end    
    if isempty(i)
        i = kMax;
    end
    %
    S.controller.ref.x = S.path.coordinates(1,S.path.acumIndx);
    S.controller.ref.y = S.path.coordinates(2,S.path.acumIndx);
    dx                 = S.path.coordinates(1,S.path.acumIndx)-S.path.coordinates(1,S.path.acumIndx-1);
    dy                 = S.path.coordinates(2,S.path.acumIndx)-S.path.coordinates(2,S.path.acumIndx-1);
    new_angle           = atan2c(dx, dy, S.controller.angle0_old);
    S.controller.ref.theta  = new_angle;
    S.controller.angle0_old = new_angle;
    S.controller.ref.xyth   = [S.controller.ref.x;S.controller.ref.y;S.controller.ref.theta];
    %
    S.path.vDir     = vDir./norm(vDir);
    vTrToObs        = S.data.xest(2*S.config.Nt+2:2*S.config.Nt+3,end);
    S.path.vTrToObs = S.path.coordinates(:,S.path.acumIndx)-vTrToObs;
    S.path.vTrToObs = S.path.vTrToObs./norm(S.path.vTrToObs);
    S.path.vDirDotvTrToObs = dot(S.path.vDir, S.path.vTrToObs);
end

function S = handleMovingTarget2(S)
    kMax        = length(S.path.coordinates)-1;
    vTrToObs    = S.data.xest(2*S.config.Nt+2:2*S.config.Nt+3,end);
    for i=S.path.acumIndx+1:kMax
        vDir                    = S.path.coordinates(:,i)-S.path.coordinates(:,i-1);
        S.path.acumIndx         = i;
        S.path.vDir             = vDir./norm(vDir);        
        S.path.vTrToObs         = S.path.coordinates(:,S.path.acumIndx)-vTrToObs;
        S.path.vTrToObs         = S.path.vTrToObs./norm(S.path.vTrToObs);
        S.path.vDirDotvTrToObs  = dot(S.path.vDir, S.path.vTrToObs);

        if S.path.vDirDotvTrToObs>0
            break;
        end        
    end
    S.controller.ref.x = S.path.coordinates(1,S.path.acumIndx);
    S.controller.ref.y = S.path.coordinates(2,S.path.acumIndx);
    dx0                 = S.path.coordinates(1,S.path.acumIndx)-S.path.coordinates(1,S.path.acumIndx-1);
    dy0                 = S.path.coordinates(2,S.path.acumIndx)-S.path.coordinates(2,S.path.acumIndx-1);
    new_angle           = atan2c(dx0, dy0, S.controller.angle0_old);
    S.controller.ref.theta  = new_angle;
    S.controller.angle0_old = new_angle;
    S.controller.ref.xyth   = [S.controller.ref.x;S.controller.ref.y;S.controller.ref.theta];
    %
    for i=1:size(S.path.listOfObs,1)
        r = (S.controller.ref.x-S.path.listOfObs(i,1))^2+(S.controller.ref.y-S.path.listOfObs(i,2))^2;
        if r <= S.path.listOfObs(i,3)^2
            [x,y]              = findOptionalPath(S,S.path.listOfObs(i,1),S.path.listOfObs(i,2),S.controller.ref.x,S.controller.ref.y,S.path.listOfObs(i,3),S.path.listOfObs(i,3));
            S.controller.ref.x = x;
            S.controller.ref.y = y;
            %
            S.controller.modVelFac = 0.25;
            flgObs = true;
            %
            break;
        end
    end
end

function S = handleTargetVelocity(S)
    flgObs = false;
    for i=1:size(S.path.listOfObs,1)
        r = (S.controller.ref.x-S.path.listOfObs(i,1))^2+(S.controller.ref.y-S.path.listOfObs(i,2))^2;
        if r <= S.path.listOfObs(i,3)^2
            [x,y]              = findOptionalPath(S,S.path.listOfObs(i,1),S.path.listOfObs(i,2),S.controller.ref.x,S.controller.ref.y,S.path.listOfObs(i,3),S.path.listOfObs(i,3));
            S.controller.ref.x = x;
            S.controller.ref.y = y;
            %
            S.controller.modVelFac = 0.25;
            flgObs = true;
            %
            break;
        end
    end
    if ~flgObs
        S.controller.modVelFac = (0.85*exp(-0.01*S.mpc.mpcCasadi.Jnum)+0.15);
    end
    if S.path.vDirDotvTrToObs < 0.1
        S = handleMovingTarget2(S);
    end
end

function S = handleWeightinMatrices(S)
     % Update matrices entiries
    R = S.mpc.Mtxs.R.*exp(-0.01*S.mpc.mpcCasadi.Jnum);
    setMtxR(S.mpc.mpcCasadi,R);
%     Q = S.mpc.Mtxs.Q./S.controller.modVelFac;
%     setMtxQ(S.mpc.mpcCasadi,Q);
    QN = S.mpc.Mtxs.QN./S.controller.modVelFac;
    setMtxQN(S.mpc.mpcCasadi,QN);
end

function S = simulateMovingObstacles(S)
% It takes the list of obstacles and changes to coordinates of one of them
    if ~isempty(S.path.listOfObs)
        for i=1:S.config.maxMovingObs
            switch i
                case 1
                    S.path.listOfObs(S.config.maxStaticObs+i,1:2) = [3+4*cos(pi+0.0675*sum(S.config.Ts*S.config.iters)), 5+3*sin(pi+0.0675*sum(S.config.Ts*S.config.iters))]; % LOS INDICES DE ESTOS OBSTACLUS DEBEN SER MENORES O IGUAL QUE EL MAX NUMBER OF STATIC OBSTACLE
                case 2
                    a = 4;
                    c = 8;
                    b = 1;
                    if S.config.accelOneMovObs
                        pot = 1.2;
                    else
                        pot = 1;
                    end
                    t = 0.015*sum(S.config.Ts*S.config.iters)^(pot);
                    x = (a*sqrt(2).*cos(-t+pi/4))./(sin(-t+pi/4).^2+1)+5.5;
                    y = (c*sqrt(2).*cos(-t+pi/4).*sin(-t+pi/4))./(sin(-t+pi/4).^2 + b)+4.5;
                    S.path.listOfObs(S.config.maxStaticObs+i,1:2) = [x, y];   % LOS INDICES DE ESTOS OBSTACLUS DEBEN SER MENORES O IGUAL QUE EL MAX NUMBER OF STATIC OBSTACLE
                case 3
                    S.path.listOfObs(S.config.maxStaticObs+i,1) = 5+6*cos(pi+0.0675*sum(S.config.Ts*S.config.iters));
                case 4
                    S.path.listOfObs(S.config.maxStaticObs+i,2) = 5+3*sin(pi+0.0675*sum(S.config.Ts*S.config.iters));
                case 5
                    S.path.listOfObs(S.config.maxStaticObs+i,1:2) = [9+1.5*cos(pi+0.08*sum(S.config.Ts*S.config.iters)), 4+1.5*sin(pi+0.08*sum(S.config.Ts*S.config.iters))];
                case 6
                    S.path.listOfObs(S.config.maxStaticObs+i,1:2) = [2+1.5*cos(pi-0.08*sum(S.config.Ts*S.config.iters)), 6+1.5*sin(pi-0.08*sum(S.config.Ts*S.config.iters))];
                case 7

                case 8

                case 9

                case 10

                case 11

                case 12
                   
                case 13

            end            
        end               
    end
end

function S = updateNearObs(S)
    if isempty(S.path.staticObs)  && isempty(S.path.dynamicObs)
        setXYRobstacle(S.mpc.mpcCasadi, [1e6, 1e6, 1]); % no obstacle found
    else
        q_k = S.data.xest(:,end);
        obsAux = [];
        for i=1:S.config.Nt+1
            pos_i               = q_k(2*S.config.Nt+1+(i-1)*2+1:2*S.config.Nt+1+i*2);
%             dis_i               = pdist2(S.path.staticObs(:,1:2),pos_i');
            dis_i               = pdist2(S.path.listOfObs(:,1:2),pos_i');
            indxObsClose        = find(dis_i<=S.config.maxDistStatObs);
            S.path.nearObs{i}   = S.path.listOfObs(indxObsClose,:);
            if ~isempty(S.path.nearObs{i})
%                 obsAux = [obsAux; S.path.nearObs{i}(1:min(size(S.path.nearObs{i},1),S.config.numStObsSgmnt),:)];
                obsAux = [obsAux; S.path.nearObs{i}(1:min(size(S.path.nearObs{i},1),S.config.totNumObs),:)];
            end
        end
        if ~isempty(obsAux)
            [B, ~, ib]      = unique(obsAux, 'rows');
            indices         = accumarray(ib, find(ib), [], @(rows){rows});
            obsAux2 = [];
            for i=1:length(indices)
                obsAux2 = [obsAux2;obsAux(indices{i}(1),:)];
            end
            %
            S.path.obsToMPC = obsAux2;
            nroObs          = size(S.path.obsToMPC,1);
            setXYRobstacle(S.mpc.mpcCasadi, S.path.obsToMPC);
            %
        end        
    end
    %
% S.path.dynObsAux = S.path.staticObs; % VER BIEN DONDE COLOAR ESTA ASIGNACION CUANDO DETECTE OBSTACULOS CON EL LIDAR
end

function S = handleStaticAndDynamicObs(S)
%     S.path.staticObs = S.path.listOfObsStr{end};
    %
    if numel(S.path.listOfObsStr)<S.config.nroFrDetMObs        
        return;
    end
    % *********************************************************************
    % Here I compute how much the obstacles have moved in succesive
    % frames.If any reachs some treshould, it could be a moving obstacle
    diagDists = [];
    for i=S.config.nroFrDetMObs-1:-1:1
        dists       = pdist2(S.path.listOfObsStr{end-i}(:,1:2), S.path.listOfObsStr{end-i+1}(:,1:2));
        diagAux     = diag(dists);
        diagDists   = [diagDists, diagAux(1:S.config.totNumObs)];
    end
    [f1,~] = find(diagDists(:,end)>S.config.thrMovingObs); % From all detected obstacles, this variable contains the indices of those that have moved.
    % Are there more moving obstacles than those that are supposed to?
    indx = f1;    
    while numel(indx)>S.config.maxMovingObs
        [v,i] = min(diagDists(indx,end));
        indx(i) = [];
    end
    S.path.dynamicObs = S.path.listOfObsStr{end}(indx,:); % update coordinates and radii on moving obstacles
    % Once moving obstacles were identified, the remaining are classified
    % as statics
    indxAux           = 1:numel(S.path.listOfObsStr{end}(:,1));
    indxStatic        = setdiff(indxAux', indx');
    S.path.staticObs  = S.path.listOfObsStr{end}(indxStatic,:); % after isolating the dynamic obstacles, the list of static obstacles is updated
    % *********************************************************************
    % Find the ellipse's equations that model each moving obstacle
    % *********************************************************************
    if ~isempty(S.path.dynamicObs)
        if isempty(S.path.dynObsEllipsePrms.a)  % If moving obstacles were found but I do not have yet modelled the ellipses
            for i=1:size(S.path.dynamicObs,1)
                indxMovObs  = pdist2(S.path.listOfObsStr{end}(:,1:2),S.path.dynamicObs(i,1:2));
                indx        = find(indxMovObs==0);
                xy          = [];
                for j=S.config.numMeasEll-1:-1:0
                    xy = [xy; S.path.listOfObsStr{end-j}(indx,1:2)];
                end
                S                               = findCoeffEllipse(S,xy(:,1),xy(:,2));
                S.path.dynObsEllipsePrms.a      = [S.path.dynObsEllipsePrms.a; S.path.dynObsEllipsePrmsAux.a];
                S.path.dynObsEllipsePrms.b      = [S.path.dynObsEllipsePrms.b; S.path.dynObsEllipsePrmsAux.b];
                S.path.dynObsEllipsePrms.xc     = [S.path.dynObsEllipsePrms.xc; S.path.dynObsEllipsePrmsAux.xc];
                S.path.dynObsEllipsePrms.yc     = [S.path.dynObsEllipsePrms.yc; S.path.dynObsEllipsePrmsAux.yc];
                S.path.dynObsEllipsePrms.t0     = [S.path.dynObsEllipsePrms.t0; S.path.dynObsEllipsePrmsAux.t0];
                S.path.dynObsEllipsePrms.kt     = [S.path.dynObsEllipsePrms.kt; S.path.dynObsEllipsePrmsAux.kt];
                S.path.dynObsEllipsePrms.theta  = [S.path.dynObsEllipsePrms.theta; S.path.dynObsEllipsePrmsAux.theta];
                S.path.dynObsEllipsePrms.t      = [S.path.dynObsEllipsePrms.t; S.path.dynObsEllipsePrmsAux.t];
                S.path.dynObsEllipsePrms.tOld   = [S.path.dynObsEllipsePrms.tOld; S.path.dynObsEllipsePrmsAux.tOld];
            end            
            return;
        elseif size(S.path.dynObsEllipsePrms.a,1)>size(S.path.dynamicObs,1) % If I currently have more ellipses than moving obstacles
            while size(S.path.dynObsEllipsePrms.a,1)>size(S.path.dynamicObs,1)
                D = [];
                for i=1:size(S.path.dynamicObs,1)
                    d = [];
                    for j=1:size(S.path.dynObsEllipsePrms.a,1)
                        xObs = S.path.dynamicObs(i,1);
                        yObs = S.path.dynamicObs(i,2);
                        d    = [d; sqrt( (xObs-S.path.dynObsEllipsePrms.xc(j)-S.path.dynObsEllipsePrms.a(j)*cos(S.path.dynObsEllipsePrms.t(i)))^2+...
                                         (yObs-S.path.dynObsEllipsePrms.yc(j)-S.path.dynObsEllipsePrms.b(j)*sin(S.path.dynObsEllipsePrms.t(i)))^2  )];
                    end
                    D = [D, d];
                end
                if size(D,2)==1
                    [~,indx] = max(D);                    
                else                    
                    [~,indx] = max(min(D));
%                     indx     = setdiff(size(D,1),indx');
%                     if numel(indx)>1
%                         indx = indx(1);
%                     end

                end
                S.path.dynObsEllipsePrms.a(indx)    = [];
                S.path.dynObsEllipsePrms.b(indx)    = [];
                S.path.dynObsEllipsePrms.xc(indx)   = [];
                S.path.dynObsEllipsePrms.yc(indx)   = [];
                S.path.dynObsEllipsePrms.t0(indx)   = [];
                S.path.dynObsEllipsePrms.kt(indx)   = [];
                S.path.dynObsEllipsePrms.theta(indx)= [];
                S.path.dynObsEllipsePrms.t(indx)    = [];
                S.path.dynObsEllipsePrms.tOld(indx) = [];
            end
        else % If I have more dynamic obstacles than ellipses
            while size(S.path.dynObsEllipsePrms.a,1)<size(S.path.dynamicObs,1)
                for i=1:size(S.path.dynamicObs,1)
                    % d = [];
                    for j=1:size(S.path.dynObsEllipsePrms.a,1) % On the go, I check if equations still a good approximation
                        xObs = S.path.dynamicObs(i,1);
                        yObs = S.path.dynamicObs(i,2);
                        d = sqrt( (xObs-S.path.dynObsEllipsePrms.xc(j)-S.path.dynObsEllipsePrms.a(j)*cos(S.path.dynObsEllipsePrms.t(j)))^2+...
                                  (yObs-S.path.dynObsEllipsePrms.yc(j)-S.path.dynObsEllipsePrms.b(j)*sin(S.path.dynObsEllipsePrms.t(j)))^2);
                        if d>S.config.distToUpdEll
                            xy   = [];
                            for k=S.config.numMeasEll-1:-1:0
                                xy = [xy; S.path.listOfObsStr{end-k}(i,1:2)];
                            end
                            S = findCoeffEllipse(S,xy(:,1),xy(:,2));
                            S.path.dynObsEllipsePrms.a(j)      = S.path.dynObsEllipsePrmsAux.a;
                            S.path.dynObsEllipsePrms.b(j)      = S.path.dynObsEllipsePrmsAux.b;
                            S.path.dynObsEllipsePrms.xc(j)     = S.path.dynObsEllipsePrmsAux.xc;
                            S.path.dynObsEllipsePrms.yc(j)     = S.path.dynObsEllipsePrmsAux.yc;
                            S.path.dynObsEllipsePrms.t0((j-1)*S.config.numMeasEll+1:j*S.config.numMeasEll) = S.path.dynObsEllipsePrmsAux.t0;
                            S.path.dynObsEllipsePrms.kt(j)     = S.path.dynObsEllipsePrmsAux.kt;
                            S.path.dynObsEllipsePrms.theta(j)  = S.path.dynObsEllipsePrmsAux.theta;
                            S.path.dynObsEllipsePrms.t(j)      = S.path.dynObsEllipsePrmsAux.t;
                            S.path.dynObsEllipsePrms.tOld(j)   = S.path.dynObsEllipsePrmsAux.tOld;
                        end
                    end
                end % May happen that all ellipses describes well the new moving obstacles added to the list. In this case. the former way of trying add a new ellipse may fail.
                % if size(S.path.dynObsEllipsePrms.a,1)<size(S.path.dynamicObs,1)
                xy   = []; % Here is where I effectively add tthe missing equation
                for k=S.config.numMeasEll-1:-1:0
                    xy = [xy; S.path.listOfObsStr{end-k}(i,1:2)];
                end
                S = findCoeffEllipse(S,xy(:,1),xy(:,2));
                S.path.dynObsEllipsePrms.a      = [S.path.dynObsEllipsePrms.a; S.path.dynObsEllipsePrmsAux.a];
                S.path.dynObsEllipsePrms.b      = [S.path.dynObsEllipsePrms.b; S.path.dynObsEllipsePrmsAux.b];
                S.path.dynObsEllipsePrms.xc     = [S.path.dynObsEllipsePrms.xc; S.path.dynObsEllipsePrmsAux.xc];
                S.path.dynObsEllipsePrms.yc     = [S.path.dynObsEllipsePrms.yc; S.path.dynObsEllipsePrmsAux.yc];
                S.path.dynObsEllipsePrms.t0     = [S.path.dynObsEllipsePrms.t0; S.path.dynObsEllipsePrmsAux.t0];
                S.path.dynObsEllipsePrms.kt     = [S.path.dynObsEllipsePrms.kt; S.path.dynObsEllipsePrmsAux.kt];
                S.path.dynObsEllipsePrms.theta  = [S.path.dynObsEllipsePrms.theta; S.path.dynObsEllipsePrmsAux.theta];
                S.path.dynObsEllipsePrms.t      = [S.path.dynObsEllipsePrms.t; S.path.dynObsEllipsePrmsAux.t];
                S.path.dynObsEllipsePrms.tOld   = [S.path.dynObsEllipsePrms.tOld; S.path.dynObsEllipsePrmsAux.tOld];
                % end
            end
        end
        % Now I have the same number of dynamic obstacles and ellipses
        indxObs = 0;   % I have to check if the current ellipses are good approximation to the moving obstacles, otherwise, I have to update the ellipse's equations
        D       = [];
        t0Temp  = [];
        for j=1:size(S.path.dynObsEllipsePrms.a,1)
            d = [];
            for i=1:size(S.path.dynamicObs,1)
                xObs = S.path.dynamicObs(i,1);
                yObs = S.path.dynamicObs(i,2);
                dis = sqrt( (xObs-S.path.dynObsEllipsePrms.xc(j)-S.path.dynObsEllipsePrms.a(j)*cos(S.path.dynObsEllipsePrms.t(j)))^2+...
                            (yObs-S.path.dynObsEllipsePrms.yc(j)-S.path.dynObsEllipsePrms.b(j)*sin(S.path.dynObsEllipsePrms.t(j)))^2  );
                if any(dis<S.config.distToUpdEll)
                    d = [d;0];
                    t0Temp = [];
                else
                    d = [d;dis];
                end                
            end
            D = [D, d];
        end
        [f,~] = find(D==0);
        f     = unique(f);
        if isempty(f)
            S.path.dynObsEllipsePrms.a      = [];
            S.path.dynObsEllipsePrms.b      = [];
            S.path.dynObsEllipsePrms.xc     = [];
            S.path.dynObsEllipsePrms.yc     = [];
            S.path.dynObsEllipsePrms.t0     = [];
            S.path.dynObsEllipsePrms.kt     = [];
            S.path.dynObsEllipsePrms.theta  = [];
            S.path.dynObsEllipsePrms.t      = [];
            S.path.dynObsEllipsePrms.tOld   = [];
            for i=1:size(S.path.dynamicObs,1)
                indxMovObs  = pdist2(S.path.listOfObsStr{end}(:,1:2),S.path.dynamicObs(i,1:2));
                indx        = find(indxMovObs==0);
                xy          = [];
%                 indxsEllip  = randperm(length(S.path.listOfObsStr),S.config.numMeasEll);
                for j=S.config.numMeasEll-1:-1:0
                    xy = [xy; S.path.listOfObsStr{end-(j)}(indx,1:2)];
                end
                S                               = findCoeffEllipse(S,xy(:,1),xy(:,2));
                S.path.dynObsEllipsePrms.a      = [S.path.dynObsEllipsePrms.a; S.path.dynObsEllipsePrmsAux.a];
                S.path.dynObsEllipsePrms.b      = [S.path.dynObsEllipsePrms.b; S.path.dynObsEllipsePrmsAux.b];
                S.path.dynObsEllipsePrms.xc     = [S.path.dynObsEllipsePrms.xc; S.path.dynObsEllipsePrmsAux.xc];
                S.path.dynObsEllipsePrms.yc     = [S.path.dynObsEllipsePrms.yc; S.path.dynObsEllipsePrmsAux.yc];
                S.path.dynObsEllipsePrms.t0     = [S.path.dynObsEllipsePrms.t0; S.path.dynObsEllipsePrmsAux.t0];
                S.path.dynObsEllipsePrms.kt     = [S.path.dynObsEllipsePrms.kt; S.path.dynObsEllipsePrmsAux.kt];
                S.path.dynObsEllipsePrms.theta  = [S.path.dynObsEllipsePrms.theta; S.path.dynObsEllipsePrmsAux.theta];
                S.path.dynObsEllipsePrms.t      = [S.path.dynObsEllipsePrms.t; S.path.dynObsEllipsePrmsAux.t];
                S.path.dynObsEllipsePrms.tOld   = [S.path.dynObsEllipsePrms.tOld; S.path.dynObsEllipsePrmsAux.tOld];
            end
        else
            indxTot     = [1:size(S.path.dynObsEllipsePrms.a,1)]';
            indxToLoop  = setdiff(indxTot, f');
% indxToLoop(find(indxToLoop>size(S.path.dynamicObs,1)))=[];            
            if ~isempty(indxToLoop)
                for i=indxToLoop
                    indxMovObs  = pdist2(S.path.listOfObsStr{end}(:,1:2),S.path.dynamicObs(i,1:2));
                    indx        = find(indxMovObs(:,1)==0);
                    xy          = [];
                    for j=S.config.numMeasEll-1:-1:0
                        xy = [xy; S.path.listOfObsStr{end-j}(indx,1:2)];
                    end
                    S                                  = findCoeffEllipse(S,xy(:,1),xy(:,2));
                    S.path.dynObsEllipsePrms.a(i)      = S.path.dynObsEllipsePrmsAux.a;
                    S.path.dynObsEllipsePrms.b(i)      = S.path.dynObsEllipsePrmsAux.b;
                    S.path.dynObsEllipsePrms.xc(i)     = S.path.dynObsEllipsePrmsAux.xc;
                    S.path.dynObsEllipsePrms.yc(i)     = S.path.dynObsEllipsePrmsAux.yc;
    %                     S.path.dynObsEllipsePrms.t0(i)     = S.path.dynObsEllipsePrmsAux.t0;
                    S.path.dynObsEllipsePrms.kt(i)     = S.path.dynObsEllipsePrmsAux.kt;
                    S.path.dynObsEllipsePrms.theta(i)  = S.path.dynObsEllipsePrmsAux.theta;
                    S.path.dynObsEllipsePrms.t(i)      = S.path.dynObsEllipsePrmsAux.t;
                    S.path.dynObsEllipsePrms.tOld(i)   = S.path.dynObsEllipsePrmsAux.tOld;
                end
            end
        end
        % Update the phsae for eahc moving obstacle as it moves along the
        % ellipse
        for j=1:length(S.path.dynObsEllipsePrms.a)
            xObs = S.path.dynamicObs(j,1);
            yObs = S.path.dynamicObs(j,2);
            S.path.dynObsEllipsePrms.t(j)     = atan2c((xObs-S.path.dynObsEllipsePrms.xc(j))/S.path.dynObsEllipsePrms.a(j), (yObs-S.path.dynObsEllipsePrms.yc(j))/S.path.dynObsEllipsePrms.b(j), S.path.dynObsEllipsePrms.t(j));
            S.path.dynObsEllipsePrms.tOld(j)  = S.path.dynObsEllipsePrms.t(j);                    
        end
        setT0MovingObs(S.mpc.mpcCasadi,S.path.dynObsEllipsePrms.t);
        






    elseif ~isempty(S.path.dynObsEllipsePrms.a) % If I do not have any moving obstacle
        for i=1:size(S.path.dynObsEllipsePrms.a,1)
            S.path.dynObsEllipsePrms.a      = [];
            S.path.dynObsEllipsePrms.b      = [];
            S.path.dynObsEllipsePrms.xc     = [];
            S.path.dynObsEllipsePrms.yc     = [];
            S.path.dynObsEllipsePrms.t0     = [];
            S.path.dynObsEllipsePrms.kt     = [];
            S.path.dynObsEllipsePrms.theta  = [];
            S.path.dynObsEllipsePrms.t      = [];
            S.path.dynObsEllipsePrms.tOld   = 0;
        end
    end



%     if ~isempty(S.path.dynObsAux)
% %         [~,indxs] = setdiff(S.path.listOfObs, S.path.dynObsAux,'rows');
% %         if ~isempty(indxs)
% %             S.path.dynObsXYcoords = [S.path.dynObsXYcoords; S.path.dynamicObs(indxs,1:2)];
% S.path.dynObsXYcoords = [S.path.dynObsXYcoords; S.path.dynamicObs(1,1:2)];
% %         end
%     end
%     if size(S.path.dynObsXYcoords,1) >= S.config.numMeasEll
%         d = sqrt((S.path.dynObsXYcoords(end,1)-S.path.dynObsEllipsePrms.xc-S.path.dynObsEllipsePrms.a*cos(S.path.dynObsEllipsePrms.t))^2+(S.path.dynObsXYcoords(end,2)-S.path.dynObsEllipsePrms.yc-S.path.dynObsEllipsePrms.b*sin(S.path.dynObsEllipsePrms.t))^2);
%         if any(d>0.2) || ~any(d) % When the measurement is beyond that 0.2 (m) with respect to the equation, it is computed again
% %             indxsEllip  = randperm(length(S.path.dynObsXYcoords),S.config.numMeasEll);
%             S               = findCoeffEllipse(S,S.path.dynObsXYcoords(end-S.config.numMeasEll+1:end,1),S.path.dynObsXYcoords(end-S.config.numMeasEll+1:end,2));
            setDynObsPrms(S.mpc.mpcCasadi,S.mpc.Mtxs.ampGaussDyn*ones(length(S.path.dynObsEllipsePrms.xc),1),S.path.dynObsEllipsePrms.xc,S.path.dynObsEllipsePrms.yc,S.path.dynObsEllipsePrms.kt,S.path.dynObsEllipsePrms.a,S.path.dynObsEllipsePrms.b,S.path.dynObsEllipsePrms.theta);
            %
%         end
                
%     end
end

function S = handleObstacles(S)
    if S.config.maxMovingObs > 0
        S = simulateMovingObstacles(S);
    end
    %
%     if S.config.maxStaticObs > 0
    if S.config.totNumObs > 0
        S.path.listOfObsStr{end+1} = S.path.listOfObs;
        if length(S.path.listOfObsStr)>=S.config.nroFrDetMObs
            S = handleStaticAndDynamicObs(S);
        end
        % From here, call the other functions
        S = updateNearObs(S);
    end
    % Save osbatcle positions for plotting purposes
    S.path.posObs{S.config.iters} = S.path.listOfObs;
end

function S = handleOcclusion(S)
    S.path.occluded.flg     = false;
    S.controller.modVelFac  = 1;
    q_k                     = S.data.xest(:,end);
    xyTgt                   = [S.controller.ref.x;S.controller.ref.y];
    vDir                    = xyTgt-q_k(2*S.config.Nt+2:2*S.config.Nt+3);
    vPerp                   = [-vDir(2);vDir(1)]./norm(vDir);
    if ~isempty(S.path.nearObs)
        nearObs = [S.path.nearObs{1};S.path.dynamicObs];
        for i=1:size(nearObs,1)
            [S,Jdis] = isOccluded(S,q_k(2*S.config.Nt+2),q_k(2*S.config.Nt+3),nearObs(i,1),nearObs(i,2),nearObs(i,3),vDir(1),vDir(2));
            if Jdis<1e-3
                S.path.occluded.p1      = q_k(2*S.config.Nt+2:2*S.config.Nt+3)';
                S.path.occluded.p2      = (xyTgt+vPerp)';
                S.path.occluded.p3      = (xyTgt-vPerp)';                
                S.path.occluded.flg     = true;
                S.controller.modVelFac  = 0.25;
                break;
            end
        end
    end    
    %
    %
%     S = findReachableTgt(S);
end

function S = call_PFA(S)
    tic;
    % Handle moving target ================================================
    S = updateTarget(S);
    % Handle target's velocity ============================================
    S = handleTargetVelocity(S);
    % Handle weighting matrices ===========================================
    S = handleWeightinMatrices(S);
    % Handle obstacles ====================================================
    S = handleObstacles(S);
    % Handle occlusion ====================================================
    S = handleOcclusion(S);
    % Update Omega_O & X_N_Omega_I ========================================
    S = updateSets(S);
    % Compute reachable target ============================================
    S = cmptReachableTgt(S);
    %
    S.exec_time.t_pfa = [S.exec_time.t_pfa, toc];
end

function S = cmptReachableTgt(S)
   isIn = S.config.X_N_Omega_I_updt.contains([S.controller.ref.x; S.controller.ref.y]);
   if ~isIn% || S.path.occluded.flg
       S = findReachableTgt(S);       
   else
       S.controller.ref.xReachable = S.controller.ref.x;
       S.controller.ref.yReachable = S.controller.ref.y;
       S.controller.ref.thetaReachable = S.controller.ref.theta;
   end
end

function S = makeSolverFindReachablePoint(S)
    numIneqs = 4;
    dimIneqs = 2;
%     x        = casadi.MX.sym('x');
%     y        = casadi.MX.sym('y');
%     x0       = casadi.MX.sym('x0');
%     y0       = casadi.MX.sym('y0');
%     H        = casadi.MX.sym('H',numIneqs,dimIneqs+1);
% 
%     optVar   = {x, y};
%     optPrm   = {x0, y0, vec(H)};

%     J        = (x-x0)^2 + (y-y0)^2;

%     ineqConstraints                           = {};
%     S.path.reachablePoint.ineqConstraints_lb  = -casadi.DM.inf(numIneqs);
%     S.path.reachablePoint.ineqConstraints_ub  = casadi.DM.zeros(numIneqs);
% 
%     for i=1:numIneqs
%        ineqConstraints = [ineqConstraints, {H(i,1)*x + H(i,2)*y - H(i,3)}]; 
%     end
    
%     problem     = struct('f',J,'x', vertcat(optVar{:}),'g', vertcat(ineqConstraints{:}),'p', vertcat(optPrm{:}));
%     nlpoptions                                  = struct;
%     nlpoptions.ipopt.max_iter                   = 2000;  %2000
%     nlpoptions.ipopt.print_level                = 0;
%     nlpoptions.print_time                       = 0;
%     nlpoptions.ipopt.acceptable_tol             = 1e-8;
%     nlpoptions.ipopt.acceptable_obj_change_tol  = 1e-6;
%     %
%     S.path.reachablePoint.solver = casadi.nlpsol('solver', 'ipopt', problem, nlpoptions);

    S.path.solverFindReachPoint = casadi.Opti();
    S.path.xr  = S.path.solverFindReachPoint.variable();
    S.path.yr  = S.path.solverFindReachPoint.variable();
    S.path.xi  = S.path.solverFindReachPoint.parameter(S.config.Nt+1);
    S.path.yi  = S.path.solverFindReachPoint.parameter(S.config.Nt+1);
    S.path.xj  = S.path.solverFindReachPoint.parameter((S.config.Nt+1)*(S.config.totNumObs)); % foer each segment, leave out every obstacle
    S.path.yj  = S.path.solverFindReachPoint.parameter((S.config.Nt+1)*(S.config.totNumObs));
    S.path.rj  = S.path.solverFindReachPoint.parameter((S.config.Nt+1)*(S.config.totNumObs));
    S.path.x0  = S.path.solverFindReachPoint.parameter();                      % x coordinate of the original reference
    S.path.y0  = S.path.solverFindReachPoint.parameter();                      % y coordinate of the original reference
    S.path.H   = S.path.solverFindReachPoint.parameter(numIneqs,dimIneqs+1);   % dimIneqs+1=3 because it is intended for a region in the 2D plane
    %
    S.path.solverFindReachPoint.set_value(S.path.H,S.config.X_N_Omega_I_updt.H);
    S.path.solverFindReachPoint.set_value(S.path.xj,1e6.*ones((S.config.Nt+1)*(S.config.totNumObs),1));
    S.path.solverFindReachPoint.set_value(S.path.yj,1e6.*ones((S.config.Nt+1)*(S.config.totNumObs),1));
    S.path.solverFindReachPoint.set_value(S.path.rj,ones((S.config.Nt+1)*(S.config.totNumObs),1));
    %
    J        = (S.path.xr-S.path.x0)^2 + (S.path.yr-S.path.y0)^2;
    %
    S.path.solverFindReachPoint.subject_to( S.path.H(:,1).*S.path.xr + S.path.H(:,2).*S.path.yr <= S.path.H(:,3) )
    %
    for i=1:S.config.Nt+1
        for j=1:S.config.totNumObs
            for k=0:1/5:1
                xAux = S.path.xi(i)+(S.path.xr-S.path.xi(i))*k;
                yAux = S.path.yi(i)+(S.path.yr-S.path.yi(i))*k;
                S.path.solverFindReachPoint.subject_to( S.path.rj(j)^2 < (xAux-S.path.xj(j))^2 + (yAux-S.path.yj(j))^2 < inf )
            end
        end
    end


    S.path.solverFindReachPoint.minimize(J);
    %
    p_optsMPC = struct('expand',true);
    s_optsMPC = struct('sb','yes','print_level',0,'gamma_theta',1e-1);%,'jacobian_approximation','exact','fast_step_computation','yes','warm_start_init_point','yes'); % 
    S.path.solverFindReachPoint.solver('ipopt',p_optsMPC,s_optsMPC);
    %
    %            
end

function S = findReachableTgt(S)
%     sol  = S.path.reachablePoint.solver('x0',[],'lbx',-casadi.DM.inf(2),'ubx',casadi.DM.inf(2),'lbg',S.path.reachablePoint.ineqConstraints_lb,'ubg',S.path.reachablePoint.ineqConstraints_ub,...
%         'p',vertcat(S.controller.ref.x,S.controller.ref.y,vec(S.config.X_N_Omega_I_updt.H)));
%     x_opt = full(sol.x);
%     S.controller.ref.xReachable = x_opt(1);
%     S.controller.ref.yReachable = x_opt(2);
%     %
    qk = S.data.xest(:,end);
    xk = qk(2*S.config.Nt+2);
    yk = qk(2*S.config.Nt+3);    
%     S.controller.ref.thetaReachable = atan2((S.controller.ref.yReachable-yk),(S.controller.ref.xReachable-xk));


    qk = S.data.xest(:,end);
    S.path.solverFindReachPoint.set_value(S.path.xi,qk(2*S.config.Nt+2:2:end-1));
    S.path.solverFindReachPoint.set_value(S.path.yi,qk(2*S.config.Nt+3:2:end));
    xjAux = [];
    yjAux = [];
    rjAux = [];
    if isempty(S.path.nearObs)
        return
    end
    for i=1:S.config.Nt+1
        xjAux = [xjAux; S.path.nearObs{i}(:,1) ];
        yjAux = [yjAux; S.path.nearObs{i}(:,2) ];
        rjAux = [rjAux; S.path.nearObs{i}(:,3) ];
    end
    numelements = size(rjAux,1);
    S.path.solverFindReachPoint.set_value(S.path.xj,1e6.*ones((S.config.Nt+1)*(S.config.totNumObs),1));
    S.path.solverFindReachPoint.set_value(S.path.yj,1e6.*ones((S.config.Nt+1)*(S.config.totNumObs),1));
    S.path.solverFindReachPoint.set_value(S.path.rj,ones((S.config.Nt+1)*(S.config.totNumObs),1));
    S.path.solverFindReachPoint.set_value(S.path.xj(1:numelements),xjAux);
    S.path.solverFindReachPoint.set_value(S.path.yj(1:numelements),yjAux);
    S.path.solverFindReachPoint.set_value(S.path.rj(1:numelements),rjAux);
    S.path.solverFindReachPoint.set_value(S.path.x0,S.controller.ref.x);
    S.path.solverFindReachPoint.set_value(S.path.y0,S.controller.ref.y);
    S.path.solverFindReachPoint.set_value(S.path.H,S.config.X_N_Omega_I_updt.H);

    try
        S.path.solutionRefAux = S.path.solverFindReachPoint.solve();
%         S.path.solverFindReachPoint.set_initial(S.path.solutionRefAux.value_variables());
        %
        S.controller.ref.xReachable = S.path.solutionRefAux.value(S.path.xr);
        S.controller.ref.yReachable = S.path.solutionRefAux.value(S.path.yr);
        S.controller.ref.thetaReachable = atan2c((S.controller.ref.xReachable-xk),(S.controller.ref.yReachable-yk),S.controller.ref.theta);%atan2((S.controller.ref.yReachable-yk),(S.controller.ref.xReachable-xk));
    catch
        
    end
    
    

end

function S = updateSets(S)
    qk          = S.data.xest(:,end);
    theta0      = qk(S.config.Nt+1);
    xy0         = qk(2*S.config.Nt+2:2*S.config.Nt+3);
    R           = [cos(-theta0) -sin(-theta0); sin(-theta0) cos(-theta0)];
    S.config.X_N_Omega_I_updt = S.config.X_N_Omega_I*R + xy0;
    %
    R           = [cos(-S.controller.ref.theta) -sin(-S.controller.ref.theta); sin(-S.controller.ref.theta) cos(-S.controller.ref.theta)];
    S.config.Omega_O_updt     = S.config.Omega_O*R + [S.controller.ref.x; S.controller.ref.y];
end

function S = makeSolverFindAltPath(S) % parametrs of an ellipse: xc,yc,a,b
    t       = casadi.MX.sym('t');
    %
    xc      = casadi.MX.sym('xc');
    yc      = casadi.MX.sym('yc');
    x1      = casadi.MX.sym('x1');
    y1      = casadi.MX.sym('y1');
    a       = casadi.MX.sym('a');
    b       = casadi.MX.sym('b');

    optVar      = {t{:}};
    optParam    = {xc, yc, x1, y1, a, b};


    J               = sqrt((xc-x1+a*cos(t))^2 + (yc-y1+b*sin(t))^2);
    Constraints     = [{}];
    nlp_constraints = [Constraints(:)'];

    problem         = struct('f',J,'x', t,'g', [],'p', vertcat(optParam{:}));
    nlpoptions                                  = struct;
    nlpoptions.ipopt.max_iter                   = 2000;  %2000
    nlpoptions.ipopt.print_level                = 0;
    nlpoptions.print_time                       = 0;
    nlpoptions.ipopt.acceptable_tol             = 1e-8;
    nlpoptions.ipopt.acceptable_obj_change_tol  = 1e-6;
%     nlpoptions.ipopt.max_cpu_time               = 30;

    S.path.solverAltPath = casadi.nlpsol('solver', 'ipopt', problem, nlpoptions);
end

function S = isTargetOccludedMakeSolver(S)
    x0      = casadi.MX.sym('x0');
    y0      = casadi.MX.sym('y0');
    xc      = casadi.MX.sym('xc');
    yc      = casadi.MX.sym('yc');
    R       = casadi.MX.sym('R');
    vx      = casadi.MX.sym('vx');
    vy      = casadi.MX.sym('vy');
    %
    tc      = casadi.MX.sym('tc');
    tr      = casadi.MX.sym('tr');
   
    optVar      = {tc,tr};
        

    S.path.occluded_optVar_lb   = [0;0];
    S.path.occluded_optVar_ub   = [1;1];

    optParam = {x0,y0,xc,yc,R,vx,vy};

    J = (x0+vx*tr - xc-R*cos(2*pi*tc))^2 + (y0+vy*tr - yc-R*sin(2*pi*tc))^2;

    problem     = struct('f',J,'x',vertcat(optVar{:}),'g',[],'p',vertcat(optParam{:}));
    nlpoptions                                  = struct;
    nlpoptions.ipopt.max_iter                   = 2000;  %2000
    nlpoptions.ipopt.print_level                = 0;
    nlpoptions.print_time                       = 0;
    nlpoptions.ipopt.acceptable_tol             = 1e-8;
    nlpoptions.ipopt.acceptable_obj_change_tol  = 1e-6;

    S.path.isTargertOccludeSolver = casadi.nlpsol('solver', 'ipopt', problem, nlpoptions);       
end

function [S,dis] = isOccluded(S,x0,y0,xc,yc,R,vx,vy)
    sol     = S.path.isTargertOccludeSolver('x0',S.path.occluded.x0,'lbx',S.path.occluded_optVar_lb,'ubx',S.path.occluded_optVar_ub,'lbg',[],'ubg',[],'p',[x0;y0;xc;yc;R;vx;vy]);
    x_opt   = full(sol.x);
    S.path.occluded.x0 = x_opt;
    dis     = full(sol.f);
end

function S = makeSolverFindCoeffEllipse(S)
    nroObservations = S.config.numMeasEll;
    xc      = casadi.MX.sym('xc');
    yc      = casadi.MX.sym('yc');
    phase   = casadi.MX.sym('t0',nroObservations);
    kt      = casadi.MX.sym('kt');
    a       = casadi.MX.sym('a');
    b       = casadi.MX.sym('b');
    theta   = casadi.MX.sym('theta');

    xP      = casadi.MX.sym('xP',nroObservations);
    yP      = casadi.MX.sym('yP',nroObservations);

    optVar      = {xc,yc,kt,a,b,theta,phase{:}};
        

    S.path.dynObsSolver_optVar_lb   = [0;0;-1;0;0;0;repmat(-20*pi,nroObservations,1)];
    S.path.dynObsSolver_optVar_ub   = [10;10;1;10;10;0;repmat(20*pi,nroObservations,1)];


    optParam            = {xP{:},yP{:}};
    deltatConstraints   = {};
    eqConstraints       = {};

    J = 0;
    for i=1:nroObservations
        J = J + 0.5*((xP(i) - (xc + a*cos(phase{i})))^2 + (yP(i) - (yc + b*sin(phase{i})))^2 + 10*(((xP(i)-xc)/a)^2+((yP(i)-yc)/b)^2-1)^2);
%         eqConstraints = [eqConstraints, {((xP(i)-xc)/a)^2+((yP(i)-yc)/b)^2-1}];
        if i>1
            deltatConstraints = [deltatConstraints, {phase{i}-phase{i-1}}];
        end
    end
    S.path.deltatConstraints_lb = casadi.DM(repmat(-1,nroObservations-1,1));
    S.path.deltatConstraints_ub = casadi.DM(repmat(1,nroObservations-1,1));
    S.path.eqConstraints_lb = casadi.DM.zeros(nroObservations);
    S.path.eqConstraints_ub = casadi.DM.zeros(nroObservations);
    %
    S.path.Constraints_lb = [S.path.deltatConstraints_lb];
    S.path.Constraints_ub = [S.path.deltatConstraints_ub];
    
%     Constraints = [{}];
%     nlp_constraints = [Constraints(:)'];

    problem     = struct('f',J,'x', vertcat(optVar{:}),'g', vertcat(deltatConstraints{:}),'p', vertcat(optParam{:}));
    nlpoptions                                  = struct;
    nlpoptions.ipopt.max_iter                   = 2000;  %2000
    nlpoptions.ipopt.print_level                = 0;
    nlpoptions.print_time                       = 0;
    nlpoptions.ipopt.acceptable_tol             = 1e-8;
    nlpoptions.ipopt.acceptable_obj_change_tol  = 1e-6;

    S.path.dynObsSolver = casadi.nlpsol('solver', 'ipopt', problem, nlpoptions);    

%     qpoptions = struct;
%     qpoptions.error_on_fail = false;
%     S.path.dynObsSolver = casadi.qpsol('solver', 'qpoases', problem, qpoptions);
end

function S = makeSolverFindLipschitzNcConstant(S)    
% First opt problem, which computes trajectories for then computing the
% Lypschitz constant
    q0                              = casadi.MX.sym('q0',S.system.nq);
    Qk1                             = {};
    Uk1                             = {};
    Qk2                             = {};
    Uk2                             = {};
    for i=1:S.config.Nc
        Qk1  = [Qk1(:)', {casadi.MX.sym(['Q1_{k+' num2str(i),'}'], S.system.nq)}];
        Uk1  = [Uk1(:)', {casadi.MX.sym(['U1_{k+' num2str(i-1),'}'], S.system.nu)}];
        Qk2  = [Qk2(:)', {casadi.MX.sym(['Q2_{k+' num2str(i),'}'], S.system.nq)}];
        Uk2  = [Uk2(:)', {casadi.MX.sym(['U2_{k+' num2str(i-1),'}'], S.system.nu)}];
    end        
    % Contraints of the optimisation problem
    state1Constraints                = {};
    state1_constraints_lb            = casadi.DM.zeros(S.config.Nc*S.system.nq);
    state1_constraints_ub            = casadi.DM.zeros(S.config.Nc*S.system.nq);
    
    state2Constraints                = {};
    state2_constraints_lb            = casadi.DM.zeros(S.config.Nc*S.system.nq);
    state2_constraints_ub            = casadi.DM.zeros(S.config.Nc*S.system.nq);
    %
    deltaU1Constraints               = {Uk1{1}-casadi.DM.zeros(S.system.nu)};
    deltaU1Constraints_lb            = repmat(S.mpc.box_constraints.dUluBounds(:,1),S.config.Nc,1);
    deltaU1Constraints_ub            = repmat(S.mpc.box_constraints.dUluBounds(:,2),S.config.Nc,1);

    deltaU2Constraints               = {Uk2{1}-casadi.DM.zeros(S.system.nu)};
    deltaU2Constraints_lb            = repmat(S.mpc.box_constraints.dUluBounds(:,1),S.config.Nc,1);
    deltaU2Constraints_ub            = repmat(S.mpc.box_constraints.dUluBounds(:,2),S.config.Nc,1);

    geometrical_constraints          = {};
    geometrical_constraints_lb       = zeros(S.config.Nt*2,1);
    geometrical_constraints_ub       = zeros(S.config.Nt*2,1);

    xy_acum = [0;0];
    for k=0:S.config.Nt-1
        xy_acum = xy_acum + [S.system.Lh(end-k,1)*cos(q0(S.config.Nt+k+1)) + S.system.Lh(end-k,2)*cos(q0(S.config.Nt+k+2));... 
                             S.system.Lh(end-k,1)*sin(q0(S.config.Nt+k+1)) + S.system.Lh(end-k,2)*sin(q0(S.config.Nt+k+2))];
        aux_var                     = q0(2*S.config.Nt+2:2*S.config.Nt+3) - xy_acum;
        geometrical_constraints     = [geometrical_constraints(:)', {q0(2*S.config.Nt+1+(k+1)*2+1:2*S.config.Nt+1+(k+2)*2) - aux_var }];
    end
    for i=1:S.config.Nc
        % Acommodate the constraints to shape the problem ---------
        if i==1
            Fk1                  = S.dynamic.FNt(q0,Uk1{i});
            Fk2                  = S.dynamic.FNt(q0,Uk2{i});
        else
            Fk1                  = S.dynamic.FNt(Qk1{i-1},Uk1{i});
            Fk2                  = S.dynamic.FNt(Qk2{i-1},Uk2{i});
            %
            deltaU1Constraints   = [deltaU1Constraints {Uk1{i}-Uk1{i-1}}];
            deltaU2Constraints   = [deltaU2Constraints {Uk2{i}-Uk2{i-1}}];
        end
        state1Constraints = [state1Constraints {Qk1{i}-Fk1}];
        state2Constraints = [state2Constraints {Qk2{i}-Fk2}];
        %
    end
    J = -mtimes((Qk1{end}(:)-Qk2{end}(:))', (Qk1{end}(:)-Qk2{end}(:)));
%     JBeta1 = -mtimes((Qk1{end}(1)-Qk2{end}(1))', (Qk1{end}(1)-Qk2{end}(1)));
%     Jx0 = -mtimes((Qk1{end}(2*S.config.Nt+2)-Qk2{end}(2*S.config.Nt+2))', (Qk1{end}(2*S.config.Nt+2)-Qk2{end}(2*S.config.Nt+2)));
%     Jy0 = -mtimes((Qk1{end}(2*S.config.Nt+3)-Qk2{end}(2*S.config.Nt+3))', (Qk1{end}(2*S.config.Nt+3)-Qk2{end}(2*S.config.Nt+3)));
    %
    nlp_constraints = [state1Constraints(:)',state2Constraints(:)',deltaU1Constraints(:)',deltaU2Constraints(:)',geometrical_constraints(:)'];
    nlp_constraints_lb = [state1_constraints_lb(:); state2_constraints_lb(:); deltaU1Constraints_lb(:); deltaU2Constraints_lb(:); geometrical_constraints_lb(:)];
    nlp_constraints_ub = [state1_constraints_ub(:); state2_constraints_ub(:); deltaU1Constraints_ub(:); deltaU2Constraints_ub(:); geometrical_constraints_ub(:)];
    %
    optVar = {q0,Qk1{:},Qk2{:},Uk1{:},Uk2{:}};
    optVar_lb   = vertcat( S.mpc.box_constraints.QluBounds(:,1),...
                               repmat(S.mpc.box_constraints.QluBounds(:,1), S.config.Nc, 1),...    q_k
                               repmat(S.mpc.box_constraints.QluBounds(:,1), S.config.Nc, 1),...    q_k
                               repmat(S.mpc.box_constraints.UluBounds(:,1), S.config.Nc, 1),...
                               repmat(S.mpc.box_constraints.UluBounds(:,1), S.config.Nc, 1));%       u_k
    optVar_ub   = vertcat( S.mpc.box_constraints.QluBounds(:,2),...
                               repmat(S.mpc.box_constraints.QluBounds(:,2), S.config.Nc, 1),...    q_k
                               repmat(S.mpc.box_constraints.QluBounds(:,2), S.config.Nc, 1),...    q_k
                               repmat(S.mpc.box_constraints.UluBounds(:,2), S.config.Nc, 1),...
                               repmat(S.mpc.box_constraints.UluBounds(:,2), S.config.Nc, 1));%       u_k     

    problemJ     = struct('f',J,'x',vertcat(optVar{:}),'g',vertcat(nlp_constraints{:}),'p',[]);
%     problemBeta1  = struct('f',JBeta1,'x',vertcat(optVar{:}),'g',vertcat(nlp_constraints{:}),'p',[]);
%     problemX0     = struct('f',Jx0,'x',vertcat(optVar{:}),'g',vertcat(nlp_constraints{:}),'p',[]);
%     problemY0     = struct('f',Jy0,'x',vertcat(optVar{:}),'g',vertcat(nlp_constraints{:}),'p',[]);
    nlpoptions                                  = struct;
    nlpoptions.ipopt.max_iter                   = 2000;  %2000
    nlpoptions.ipopt.print_level                = 0;
    nlpoptions.print_time                       = 0;
    nlpoptions.ipopt.acceptable_tol             = 1e-8;
    nlpoptions.ipopt.acceptable_obj_change_tol  = 1e-6;
    finddifQdifU = casadi.nlpsol('solver', 'ipopt', problemJ, nlpoptions);    
    for j=1:20
        sol = finddifQdifU('x0',S.stability.q0Opt,'lbx',optVar_lb,'ubx',optVar_ub,'lbg',nlp_constraints_lb,'ubg',nlp_constraints_ub,'p'  ,[]);        
        x_opt          = full(sol.x);            
        S.stability.q0 = reshape(x_opt(1:S.system.nq),S.system.nq,1);
        S.stability.Q1 = reshape(x_opt(S.system.nq+1:S.system.nq+S.config.Nc*S.system.nq),S.system.nq,S.config.Nc);
        S.stability.Q2 = reshape(x_opt(S.system.nq+S.config.Nc*S.system.nq+1:S.system.nq+S.config.Nc*S.system.nq+S.config.Nc*S.system.nq),S.system.nq,S.config.Nc);
        S.stability.U1 = reshape(x_opt(S.system.nq+S.config.Nc*S.system.nq+S.config.Nc*S.system.nq+1:S.system.nq+S.config.Nc*S.system.nq+S.config.Nc*S.system.nq+S.config.Nc*S.system.nu),S.system.nu,S.config.Nc);
        S.stability.U2 = reshape(x_opt(S.system.nq+S.config.Nc*S.system.nq+S.config.Nc*S.system.nq+S.config.Nc*S.system.nu+1:end),S.system.nu,S.config.Nc);
        %
        S.stability.q0Opt = x_opt;
    end

% Second opt problem, which effectively computes the Lypschitz constant
    L    = casadi.MX.sym('L');
    difQ = casadi.MX.sym('difQ',S.config.Nc);
    difU = casadi.MX.sym('difU',S.config.Nc);
    
    Lyp_constraints     = {};
    Lyp_constraints_lb  = -casadi.DM.inf(S.config.Nc);
    Lyp_constraints_ub  = casadi.DM.zeros(S.config.Nc);
    
    for i=1:S.config.Nc
        Lyp_constraints = [Lyp_constraints, {difQ(i) - L*difU(i)} ];
    end
    
    nlp_constraints     = Lyp_constraints(:)';
    nlp_constraints_lb  = [Lyp_constraints_lb(:)];
    nlp_constraints_ub  = [Lyp_constraints_ub(:)];
    
    optVar      = L;
    optVar_lb   = 0;
    optVar_ub   = inf;
    
    JBeta1           = L^2;
    optParam    = {difQ,difU};
    
    problem2     = struct('f',JBeta1,'x', optVar,'g', vertcat(nlp_constraints{:}),'p',vertcat(optParam{:}));
    nlpoptions                                  = struct;
    nlpoptions.ipopt.max_iter                   = 2000;  %2000
    nlpoptions.ipopt.print_level                = 0;
    nlpoptions.print_time                       = 0;
    nlpoptions.ipopt.acceptable_tol             = 1e-8;
    nlpoptions.ipopt.acceptable_obj_change_tol  = 1e-6;
    finddifQdifU = casadi.nlpsol('solver', 'ipopt', problem2, nlpoptions);
    %
    prmDifQ = [];
    prmDifU = [];
    for i=1:S.config.Nc
        prmDifQ = [prmDifQ; norm(S.stability.Q1(:,i)-S.stability.Q2(:,i))];
        prmDifU = [prmDifU; norm(S.stability.U1(:,i)-S.stability.U2(:,i))];
    end
    prms = vertcat(prmDifQ, prmDifU);
    %
    sol = finddifQdifU('x0',[],'lbx',optVar_lb,'ubx',optVar_ub,'lbg',nlp_constraints_lb,'ubg',nlp_constraints_ub,'p',prms);
        
    x_opt = full(sol.x);
    %
    S.stability.LipschitzConstant = x_opt;
end

function S = findCoeffEllipse(S,xObservations,yObservations)
    sol     = S.path.dynObsSolver('x0',S.path.dynObsEllipsePrmsAux.x0,'lbx',S.path.dynObsSolver_optVar_lb,'ubx',S.path.dynObsSolver_optVar_ub,'lbg',S.path.Constraints_lb,'ubg',S.path.Constraints_ub,'p',vertcat(xObservations,yObservations));
    x_opt   = full(sol.x);
    %
    S.path.dynObsEllipsePrmsAux.xc              = x_opt(1);
    S.path.dynObsEllipsePrmsAux.yc              = x_opt(2);    
%     S.path.dynObsEllipsePrmsAux.kt              = x_opt(3);
S.path.dynObsEllipsePrmsAux.kt              = x_opt(end)-x_opt(end-1);
    S.path.dynObsEllipsePrmsAux.a               = x_opt(4);
    S.path.dynObsEllipsePrmsAux.b               = x_opt(5);           
    S.path.dynObsEllipsePrmsAux.theta           = x_opt(6);
    S.path.dynObsEllipsePrmsAux.t0              = x_opt(7:end);
    %
    S.path.dynObsEllipsePrmsAux.x0              = x_opt;
    %
    S.path.dynObsEllipsePrmsAux.t = atan2c((xObservations(end)-S.path.dynObsEllipsePrmsAux.xc)/S.path.dynObsEllipsePrmsAux.a, (yObservations(end)-S.path.dynObsEllipsePrmsAux.yc)/S.path.dynObsEllipsePrmsAux.b, S.path.dynObsEllipsePrmsAux.tOld);
    S.path.dynObsEllipsePrmsAux.tOld = S.path.dynObsEllipsePrmsAux.t;
end

function [x,y] = findOptionalPath(S,xc,yc,x1,y1,a,b)
    sol     = S.path.solverAltPath('x0',0,'lbx',-2*pi,'ubx',2*pi,'lbg',0,'ubg',0,'p',vertcat(xc,yc,x1,y1,a,b));    
    x_opt   = full(sol.x);
    x       = xc + a*cos(x_opt);
    y       = yc + b*sin(x_opt);
end

function S = call_SYSTEM(S)
    if S.config.SIM == true
        if S.config.mpc
            simulation_input.x = S.data.xsim(:,end);
            simulation_input.u = S.mpc.Controls(:,end);%S.mpc.mpcCasadi.u_k;            
        else
%             simulation_input.x = [S.data.xsim(1:end-2,end); S.Michalek2017.u0];
%             simulation_input.u = [0;0];
            simulation_input.x = S.data.xsim(:,end);
            simulation_input.u = S.Michalek2017.u0;
        end
        if S.config.slip == 1
            simulation_input.slip = zeros(2*S.config.Nt+1+2*(S.config.Nt+1),1); %+7
        else
            q1                      = S.dynamic.FNt(simulation_input.x,diag(S.config.slip-[1;1])*simulation_input.u);
            q2                      = S.dynamic.FNt(simulation_input.x,[0; 0]);
            simulation_input.slip   = q1-q2;
        end
        if sum(S.config.procDist.amp) > 0
            if strcmp(S.config.procDist.type,'gaussian') || strcmp(S.config.procDist.type,'normal')
                simulation_input.dist = S.config.procDist.amp.*randn(2*S.config.Nt+1+2*(S.config.Nt+1),1); % +7
            elseif strcmp(S.config.procDist.type,'uniform')
                simulation_input.dist = S.config.procDist.amp.*(randn(2*S.config.Nt+1+2*(S.config.Nt+1),1)-1); % +7
            end        
        else
            simulation_input.dist = zeros(2*S.config.Nt+1+2*(S.config.Nt+1),1); %+7
        end
        S.data.procDist = [S.data.procDist, simulation_input.dist];
        S.data.slip     = [S.data.slip, full(simulation_input.slip)];
        %
%         states      = S.dynamic.FNt(simulation_input.x, simulation_input.u) + simulation_input.slip + simulation_input.dist;
        states      = FNt(full(simulation_input.x), full(simulation_input.u)) + simulation_input.slip + simulation_input.dist;   % generated mex file
        S.data.xsim = [S.data.xsim, full(states)];
    else
        S = write_ctrl(S, S.mpc.Husky_w0, S.mpc.Husky_v0);
    end
end

function S = call_CONTROLLER(S)
    if S.config.vNr < S.config.vNTarget
        S.config.vNr = S.config.vNr + S.config.deltavN;
    end
    if strcmp(S.config.method,'proposed')
        S = call_MPC(S);
    elseif strcmp(S.config.method,'fnmppc')
        S = call_FNMPPC(S);
    else
        error('Method not found...');
    end
end

function S = call_FNMPPC(S)
    tic;
    % Update datda
    S.fnmppc.optiMPC.set_value(S.fnmppc.qinit,S.data.xest(:,end));
    S.fnmppc.optiMPC.set_value(S.fnmppc.qref,[S.controller.ref.theta;S.controller.ref.x;S.controller.ref.y]);
    % Update obstacle information
    obsSD = [];
    if size(S.path.listOfObsStr,2)>1
        for i=1:size(S.path.listOfObsStr{end},1)
            dx      = S.path.listOfObsStr{end}(i,1)-S.path.listOfObsStr{end-1}(i,1);
            dy      = S.path.listOfObsStr{end}(i,2)-S.path.listOfObsStr{end-1}(i,2);
            obsSD   = [obsSD; [S.path.listOfObsStr{end}(i,1:3),dx,dy]];            
        end
    end
    S.fnmppc.optiMPC.set_value(S.fnmppc.ObsSD,obsSD);
    %    
    try
        S.fnmppc.solutionMPC = S.fnmppc.optiMPC.solve();
        S.fnmppc.optiMPC.set_initial(S.fnmppc.solutionMPC.value_variables());
        Utraj   = S.fnmppc.solutionMPC.value(S.fnmppc.u);
        S.fnmppc.Qtraj    = S.fnmppc.solutionMPC.value(S.fnmppc.qMPC);
    catch
        Utraj = S.fnmppc.optiMPC.debug.value(S.fnmppc.u);
        S.fnmppc.Qtraj = [];
    end    
    %    
    if S.config.iters > S.config.Ne+1        
        S.mpc.Controls  = [S.mpc.Controls, Utraj(:,1)];
    else
        S.mpc.Controls  = [S.mpc.Controls, [0;0]];
    end
    S.fnmppc.optiMPC.set_value(S.fnmppc.u_last,S.mpc.Controls(:,end));
    %
    S.mpc.Controls_w0 = [S.mpc.Controls_w0, S.mpc.Controls(1,end)];
    S.mpc.Controls_v0 = [S.mpc.Controls_v0, S.mpc.Controls(2,end)];
    % Velocities to apply to the Husky ____________________________________
    S.mpc.Husky_w0    = S.mpc.Controls_w0(end);
    S.mpc.Husky_v0    = S.mpc.Controls_v0(end);    
    %
    S.exec_time.t_mpc = [S.exec_time.t_mpc, toc];
    %
    S.data.references = [S.data.references, S.controller.ref.xyth];
end

function S = call_MPC(S)
    tic;    
    % Update data ________________________________________________________    
    setq0(S.mpc.mpcCasadi,S.data.xest(:,end));
%     setReference(S.mpc.mpcCasadi,[S.controller.ref.theta;S.controller.ref.x;S.controller.ref.y]);
    setReference(S.mpc.mpcCasadi,[S.controller.ref.thetaReachable;S.controller.ref.xReachable;S.controller.ref.yReachable]);
    % Solve control problem _______________________________________________
    solve(S.mpc.mpcCasadi);
    if S.config.iters > S.config.Ne+1        
        S.mpc.Controls  = [S.mpc.Controls, S.mpc.mpcCasadi.u_k];
    else
        S.mpc.Controls  = [S.mpc.Controls, [0;0]];
    end
    % 
    setU_l(S.mpc.mpcCasadi,S.mpc.mpcCasadi.u_k);
    S.mpc.Controls_w0 = [S.mpc.Controls_w0, S.mpc.Controls(1,end)];
    S.mpc.Controls_v0 = [S.mpc.Controls_v0, S.mpc.Controls(2,end)];
    % Velocities to apply to the Husky ____________________________________
    S.mpc.Husky_w0    = S.mpc.Controls_w0(end);
    S.mpc.Husky_v0    = S.mpc.Controls_v0(end);    
    %
    S.exec_time.t_mpc = [S.exec_time.t_mpc, toc];
    %
    S.data.references = [S.data.references, S.controller.ref.xyth];
end

function S = init_obsDetection(S)    
    %
    S.obs.pcMerge.gridSize        = 0.1;      % size of the voxels for averaging points
    S.obs.pcMerge.nroClouds       = 3;
    S.obs.pcSegdist.minNroPoints  = 15;       % minimum number of points of a cluster
    S.obs.pcSegdist.maxNroPoints  = inf;      % maximum number of points of a cluster
    S.obs.pcSegdist.minDistance   = 0.8;        % min disance between point of different clusters
    S.obs.filterWalls.area        = 2;        % area= deltaX+h. Clusters with area greater than this parameter are discarded
    S.obs.filterWalls.height      = 2;        % cluster with a height greater than this parameter are discarded
    S.obs.filterGround.disDev     = 0.2;      % deviation of points to the plane that is considered as ground
    S.obs.filterGround.angleDev   = 5;        % angular deviation of the normal vector to the plance considered as ground
    %
    S.obs.newObstacle.minDis      = 1;        % when an obstacle is found, if it is near to some other, they are merged.
    %
    S.obs.remMovingObs.searchRad  = 5;        % radius to find old obstacles
    S.obs.remMovingObs.disBetObs  = 1;        % when a new obstacle is found, and it is near to some obstacle from the list, this one is replaced by the new, assuming that it is a moving obstacle
    S.obs.remMovingObs.lossVotes  = 2;        % when osbtacles from the list are not found in the search radius, they loss votes
    %
    S.obs.validateObs.votes       = 10;       % minimum number of votes to consider a cluster as an obstacle, every time the same cluster is found, it gains a vote
    S.obs.validateObs.timeIndx    = 30;       % an obstacle with few votes and that it is not detected since a "long" time ago, is dicarded
    % variables to store pointclouds and data
    S.obs.ptCloud                 = [];
    S.obs.ptCloudAvg              = [];
    S.obs.ptCloudMap              = [];
    S.obs.Obstacles               = [];
    S.obs.ObstaclesAux            = [];
    %
    if ~exist('S.vlp16','var')
        S.vlp16 = velodynelidar('VLP16');              
    end
%     for i=1:S.obs.pcMerge.nroClouds
%         [S.obs.ptCloud, S.obs.timeStamps] = read(S.vlp16, 'latest');
%         S.obs.ptCloudAvg = [S.obs.ptCloud, S.obs.ptCloudAvg];
%     end
end

function S = call_OBSDETECTOR(S)
    tic;        
    if S.config.lidarOn
        if S.vlp16.NumPointCloudsAvailable > 0
            [S.obs.ptCloud, S.obs.timeStamps] = read(S.vlp16,'all');
            S.sensors.vlp16 = [S.sensors.vlp16, {S.obs.ptCloud}];
            if S.config.obsDetection
                S.obs.ptCloud = S.obs.ptCloud(end);
                if S.config.SIM
                    xy0     = S.data.xest(2*S.config.Nt+2:2*S.config.Nt+3,end);
                    theta0  = S.data.xest(S.config.Nt+1,end);
                else
                    xy0     = [S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0];
                    theta0  = S.ROS.sensors.vectornav_theta0;
                end
                % Remove ground
                [~,~,outlierIndicesW,~]             = pcfitplane(S.obs.ptCloud,S.obs.filterGround.disDev,[0 0 1],S.obs.filterGround.angleDev);
                S.obs.ptCloudAux                    = select(S.obs.ptCloud,outlierIndicesW);
                % Select ROI with LiDAR at centre
                S.obs.xyz                           = S.obs.ptCloudAux.Location(abs(S.obs.ptCloudAux.Location(:,1))<=S.config.lidarLimits.X & S.obs.ptCloudAux.Location(:,2)<=S.config.lidarLimits.Y & S.obs.ptCloudAux.Location(:,3)<=S.config.lidarLimits.Z,:);
                S.obs.ptCloudAux                    = pointCloud(S.obs.xyz);
                [labels,numClusters]                = pcsegdist(S.obs.ptCloudAux,S.obs.pcSegdist.minDistance,'NumClusterPoints',[S.obs.pcSegdist.minNroPoints,S.obs.pcSegdist.maxNroPoints]);
                S.obs.localObstacles                = findObstacles(S.obs.ptCloudAux, xy0, labels, numClusters, 0); % last argument passed as a constant
                S.obs.globalObstacles               = localToGlobalObs(S,S.obs.localObstacles,theta0,xy0);
                %
% S.path.listOfObs = [];    
                for i=1:size(S.obs.globalObstacles,1)    
                    S = genObstacle(S,S.obs.globalObstacles(i,1),S.obs.globalObstacles(i,2),max(S.obs.globalObstacles(i,3:4)));
                end
            end
        end       
    end
    S.exec_time.t_obsdetector = [S.exec_time.t_obsdetector, toc];    
end

function Obstacles = localToGlobalObs(S,local,theta0,xy_husky)
Obstacles = local;
    if isempty(local)
        return;
    else 
    obsOutOfLims = [];
        for i=1:size(local,1)
            d = sqrt(local(i,1)^2+local(i,2)^2); 
            alpha = atan(local(i,1)/local(i,2));    % x and y components are swapped since the lidar read in this way.
            xol = d * cos(theta0-alpha);
            yol = d * sin(theta0-alpha);
            xog = xy_husky(1) + xol;
            yog = xy_husky(2) + yol;
            if (xog >= S.config.x_min) && (xog <= S.config.x_max) && (yog >= S.config.y_min) && (yog <= S.config.y_max)
                Obstacles(i,1) = xog;
                Obstacles(i,2) = yog;
            else
                obsOutOfLims = [obsOutOfLims; i];
            end
        end   
        Obstacles(obsOutOfLims,:) = [];
    end
end

function obs = findObstacles(ptCloudAux, posH, labels, numClusters, indx)
obs         = [];
obsToDel    = [];
% ptCloudAvg = [];
    for i=1:numClusters
        obstacle_i      = find(labels == i);
        pt_obstacle     = select(ptCloudAux,obstacle_i);
% if i==1
%     ptCloudAvg = pt_obstacle;
% else
%     ptCloudAvg = pcmerge(ptCloudAvg, pt_obstacle);
% end

        x_mean          = mean(pt_obstacle.Location(:,1));
        y_mean          = mean(pt_obstacle.Location(:,2));
%         disGlo          = norm([x_mean, y_mean]);
        deltaX          = (pt_obstacle.XLimits(2)-pt_obstacle.XLimits(1));
        deltaY          = (pt_obstacle.YLimits(2)-pt_obstacle.YLimits(1));
        deltaZ          = (pt_obstacle.ZLimits(2)-pt_obstacle.ZLimits(1));        
%         h               = deltaZ;
%         area            = deltaX*h;
%         rho             = (pt_obstacle.Location(:,1)-mean(pt_obstacle.Location(:,1)))'*(pt_obstacle.Location(:,2)-mean(pt_obstacle.Location(:,2)))/length(pt_obstacle.Location)^2;% assumming here equiprobability of each point...
%         disH            = norm([x_mean-posH(1), y_mean-posH(2)]);
%         newObs          = [x_mean, y_mean, disGlo, deltaX, deltaY, rho, area, h, disH, indx, 0];
        newObs          = [x_mean, y_mean, deltaX, deltaY];
        obs             = [obs; newObs];
    end
    obsToDel        = unique(obsToDel);
    obs(obsToDel,:) = [];
% if ~isempty(ptCloudAvg)    
%     pcshow(ptCloudAvg)   
% end

end

function Obstacles = remSporaidcObs(globalObstacles, Obstacles, posH, theta, R, distanceBetObs, lossVotes)
    if isempty(Obstacles)
        return;
    end
    if ~isempty(globalObstacles)
        nroGloObs = size(globalObstacles,1);
    else
        nroGloObs = 1;
    end
    
    nroObs  = size(Obstacles,1);    
        
    for i=1:nroGloObs
        for j=1:nroObs
            posObs      = Obstacles(j,1:2)';
            if ~isempty(globalObstacles)
                posNewObs = globalObstacles(i,1:2)';
                disBetObs = sqrt((posNewObs(1)-posObs(1))^2 + (posNewObs(2)-posObs(2))^2);
            end                        
            vecObs      = posObs-posH;
            vecObs      = vecObs ./ norm(vecObs);
            dotProd     = cos(theta)*vecObs(1)+sin(theta)*vecObs(2);
            r           = sqrt((posObs(1)-posH(1))^2 + (posObs(2)-posH(2))^2);
            if r<=R && dotProd>0
                if isempty(globalObstacles)
                    Obstacles(j,11) = Obstacles(j,11) - lossVotes;
                    Obstacles(j,11) = min(Obstacles(j,11), 0);
                elseif disBetObs <= distanceBetObs
                    x_mean  = globalObstacles(i,1);
                    y_mean  = globalObstacles(i,2);
                    disGlo  = sqrt(x_mean^2+y_mean^2);
                    deltaX  = globalObstacles(i,4);
                    deltaY  = globalObstacles(i,5);
                    rho     = globalObstacles(i,6);
                    h       = globalObstacles(i,8);
                    area    = deltaX*h;
                    indx    = globalObstacles(i,10);
                    disH    = min([globalObstacles(i,9), Obstacles(j,9)]);
                    votes   = Obstacles(j,11);
                    newObs  = [x_mean, y_mean, disGlo, deltaX, deltaY, rho, area, h, disH, indx, votes];
                    Obstacles(j,:) = newObs;
                end
            end
        end
    end
end

function listObsAux = checkIfNewObs(listObsAux, recentFoundObs, minDis)
    if isempty(listObsAux)
        listObsAux = recentFoundObs;
        return;
    end
    if isempty(recentFoundObs)
        return;
    end    

    nroObsAux   = size(listObsAux,1);
    nroFoundAux = size(recentFoundObs,1);
    
    indxNewObs  = [];
    for i=1:nroFoundAux % recently found  
        countMatchs = 0;
        indxMatches = [];
        for j=1:nroObsAux % elements in the auxiliar list
            dij = sqrt((listObsAux(j,1)-recentFoundObs(i,1))^2 + (listObsAux(j,2)-recentFoundObs(i,2))^2);
            if dij <= minDis
                listObsAux(j,11) = listObsAux(j,11)+1;
                countMatchs      = countMatchs+1;
                indxMatches      = [indxMatches, [i;j]];
            else
                indxNewObs       = [indxNewObs; i];
            end
        end        
        if ~isempty(indxMatches)
            x_mean  = mean(listObsAux(indxMatches(2,:)',1));
            y_mean  = mean(listObsAux(indxMatches(2,:)',2));
            disGlo  = sqrt(x_mean^2+y_mean^2);
            deltaX  = max(listObsAux(indxMatches(2,:)',4));
            deltaY  = max(listObsAux(indxMatches(2,:)',5));
            rho     = mean(listObsAux(indxMatches(2,:)',6));
            h       = max(listObsAux(indxMatches(2,:)',8));
            area    = deltaX*h;
            indx    = max(listObsAux(indxMatches(2,:)',10));
            disH    = min(listObsAux(indxMatches(2,:)',9));
            votes   = max(listObsAux(indxMatches(2,:)',11));
            newObs  = [x_mean, y_mean, disGlo, deltaX, deltaY, rho, area, h, disH, indx, votes];
            listObsAux(indxMatches(2,:),:) = [];
            listObsAux = [listObsAux; newObs];
            nroObsAux  = size(listObsAux,1);
        end
    end
    if ~isempty(indxNewObs)
        indxNewObs = unique(indxNewObs);
        listObsAux = [listObsAux; recentFoundObs(indxNewObs,:)];  
    end
    % Sort list
    [~,indx] = sort(listObsAux(:,11));
    indx = flipud(indx);
    listObsAux = listObsAux(indx,:);
end

function StatObs = remWalls(Obstacles, areaThr, heightThr)
    StatObs  = Obstacles;
    if isempty(Obstacles)
        return;
    end    

    area     = StatObs(:,7);
    h        = StatObs(:,8);
    I        = length(area);

    obsToDel = [];
    for i=1:I
        if area(i) > areaThr || h(i) >= heightThr % then remove the farest obstacle
            obsToDel = [obsToDel; i];
        end
    end
% a= obsToDel    
    StatObs(obsToDel,:) = [];
end

function S = compute_curvature(S)
    S.path.curvature = [];
    for i=3:length(S.path.coordinates)
        x1 = S.path.coordinates(1,i);
        x2 = S.path.coordinates(1,i-1);
        x3 = S.path.coordinates(1,i-2);
        y1 = S.path.coordinates(2,i);
        y2 = S.path.coordinates(2,i-1);
        y3 = S.path.coordinates(2,i-2);
        %
        a = sqrt((x1-x2)^2 + (y1-y2)^2);
        b = sqrt((x2-x3)^2 + (y2-y3)^2);
        c = sqrt((x3-x1)^2 + (y3-y1)^2);
        s = (a+b+c)/2;
        A = sqrt(s*(s-a)*(s-b)*(s-c));
        %
        S.path.curvature = [S.path.curvature, 4*A/(a*b*c)];
    end
%     figure;plot3(S.path.coordinates(1,3:end),S.path.coordinates(2,3:end),1./S.path.curvature,'r'); grid on; zlim([0 10]);
end

function S = gen_path(S,type)
    if strcmp(type,'infinity')
        a = 4;
        c = 8;
        b = 1;
        t = 0:0.00025:2*pi;
        x = (a*sqrt(2).*cos(t))./(sin(t).^2+1);
        y = (c*sqrt(2).*cos(t).*sin(t))./(sin(t).^2 + b);
        %
        S.path.coorection_x = 5.5;  % correction for the field experiemnts in order to fit the path in my local reference frame
        S.path.coorection_y = 4.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = (a*sqrt(2).*cos(S.path.s))./(sin(S.path.s).^2+1) + S.path.coorection_x;
        S.path.fy           = (c*sqrt(2).*cos(S.path.s).*sin(S.path.s))./(sin(S.path.s).^2 + b) + S.path.coorection_y;
        %
        S = compute_curvature(S);
    elseif strcmp(type,'flat_infinity')
        a = 7;
        b = 7.5;
        t = 0:0.0005:2*pi;
        x = a*sin(t);
        y = b*sin(t).^2.*cos(t);
        %
        S.path.coorection_x = 6;  % correction for the field experiemnts in order to fit the path in my locala reference frame
        S.path.coorection_y = 4.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = a*sin(S.path.s) + S.path.coorection_x;
        S.path.fy           = b*sin(S.path.s).^2*cos(S.path.s) + S.path.coorection_y;
        %
        % fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        % fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        %
        % dfxdt               = fx_fun.jacobian;
        % dfydt               = fy_fun.jacobian;
% 
        % % Numerical approximation of the function that determines the att.
        % % val.
        % di                  = 0.01;
        % I                   = -pi:di:3*pi;
        % alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        % theta               = alpha0;
        % integrando          = [];
        % args                = [];
        % 
        % for i=I
        %     arg         = (dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])) / (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[]));
        %     args        = [args, arg];
        % %     integrando  = [integrando, atan(arg)];
        %     integrando  = [integrando, atan2((dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])), (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[])))];
        %     theta       = [theta, theta(end)+integrando(end)];
        % end               
        % % Polyfit
        % p = polyfit(I,full(theta(1:end-1)),20);
        S.path.st            = casadi.MX.sym('st');
        % S.path.ftheta        = [];%sum(S.path.s.*p);
        % for i=1:length(p)
        %     if i==1
        %         S.path.ftheta = -p(i)*S.path.st^(length(p)-i);
        %     else
        %         S.path.ftheta = S.path.ftheta - p(i)*S.path.st^(length(p)-i);
        %     end            
        % end
        %
        S.path.ftheta        = 0;%sum(S.path.s.*p);        
        S.path.r      = [ S.path.fx;  
                          S.path.fy;
                          S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type ,'segment')
        x0 = 1;
        xf = 23;
        dt = 0.01;
        t = 0:dt:2*pi;
        x = x0 + (xf-x0) * t / (2*pi);
        y0 = 6;
        y = repmat(y0,1,length(x));
        %
        S.path.coordinates  = [x ; y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.st           = casadi.MX.sym('st');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = x0 + (xf-x0) * S.path.s / (2*pi);
        S.path.fy           = y0;
        %
        fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        % Numerical approximation of the function that determines the att.
        % val.
        di                  = 0.01;
        alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        theta               = alpha0;
        S.path.ftheta       = theta;
        %
        S.path.r            = [ S.path.fx;  
                                S.path.fy;
                                S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type ,'invert_segment')
        x0 = 23;
        xf = 1;
        dt = 0.01;
        t = 0:dt:2*pi;
        x = x0 + (xf-x0) * t / (2*pi);
        y0 = 6;
        y = repmat(y0,1,length(x));
        %
        S.path.coordinates  = [x ; y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.st           = casadi.MX.sym('st');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = x0 + (xf-x0) * S.path.s / (2*pi);
        S.path.fy           = y0;
        %
        fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        % Numerical approximation of the function that determines the att.
        % val.
        di                  = 0.01;
        alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        theta               = alpha0;
        S.path.ftheta       = theta;
        %
        S.path.r            = [ S.path.fx;  
                                S.path.fy;
                                S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type,'rectangular')
        t = 0:0.01:2*pi;
        p = 3;
        q = 2.5;

        x = p.*(sqrt(cos(t).*cos(t)).*cos(t) + sqrt(sin(t).*sin(t)).*sin(t));
        y = q.*(sqrt(cos(t).*cos(t)).*cos(t) - sqrt(sin(t).*sin(t)).*sin(t));

        S.path.coorection_x = p+4;  % correction for the field experiemnts in order to fit the path in my locala reference frame
        S.path.coorection_y = q+1.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = p.*(sqrt(cos(S.path.s).*cos(S.path.s)).*cos(S.path.s) + sqrt(sin(S.path.s).*sin(S.path.s)).*sin(S.path.s)) + S.path.coorection_x;
        S.path.fy           = q.*(sqrt(cos(S.path.s).*cos(S.path.s)).*cos(S.path.s) - sqrt(sin(S.path.s).*sin(S.path.s)).*sin(S.path.s)) + S.path.coorection_y;
        %
        fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        %
        dfxdt               = fx_fun.jacobian;
        dfydt               = fy_fun.jacobian;

        % Numerical approximation of the function that determines the att.
        % val.
        di                  = 0.01;
        I                   = -pi:di:3*pi;
        alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        theta               = alpha0;
        integrando          = [];
        args                = [];
        
        for i=I
            arg         = (dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])) / (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[]));
            args        = [args, arg];
        %     integrando  = [integrando, atan(arg)];
            integrando  = [integrando, atan2((dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])), (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[])))];
            theta       = [theta, theta(end)+integrando(end)];
        end               
        % Polyfit
        p = polyfit(I,full(theta(1:end-1)),20);
        S.path.st            = casadi.MX.sym('st');
        S.path.ftheta        = [];%sum(S.path.s.*p);
        for i=1:length(p)
            if i==1
                S.path.ftheta = -p(i)*S.path.s^(length(p)-i);
            else
                S.path.ftheta = S.path.ftheta - p(i)*S.path.s^(length(p)-i);
            end            
        end
        %
        S.path.r      = [ S.path.fx;  
                          S.path.fy;
                          S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type,'test')
        % TEST TRAJECTORY #################################################
%         t = 0:3500;
%         x = 3*sin(2*pi.*t/t(end)) + 10;
%         y = 6*cos(2*pi.*t/t(end)) + 10;
%         S.path.coordinates = [x;y];
%         %
        N  = (S.config.tf-S.config.t0)/S.config.Ts;
        t           = linspace(0,2*pi,N);
        x           = 10.*cos(t)+5.*cos(5.*t) - 10;
        y           = 10.*sin(t)+5.*sin(5.*t) + 5;
        S.path.coordinates = [x;y];
    elseif strcmp(type,'monaco')
        load('monaco_xy');
        y = y.*0.7;
        S.path.coordinates  = [x(1:545); y(1:545)];
        
    elseif strcmp(type,'test_betas')
        x = 2:0.1:7.5;
        y = repmat(0.5,1,length(x));

        t = -pi/2:0.1:pi/2;
        x1 = 1.5.*cos(t) + 7.5;
        y1 = 1.5.*sin(t) + 2;

        x = [x, x1];
        y = [y, y1];

        x1 = x(end):-0.1:5;
        y1 = repmat(y(end),1,length(x1));

        x = [x, x1];
        y = [y, y1];

        t = -pi/2:-0.1:-3*pi/2;
        x1 = 1.5.*cos(t) + 4.5;
        y1 = 1.5.*sin(t) + 5;

        x = [x, x1];
        y = [y, y1];

        x1 = x(end):0.1:7.5;
        y1 = repmat(y(end),1,length(x1));

        x = [x, x1];
        y = [y, y1];

        t = -pi/2:0.1:pi/2;
        x1 = 1.5.*cos(t) + 7.5;
        y1 = 1.5.*sin(t) + 8;

        x = [x, x1];
        y = [y, y1];

        x1 = x(end):-0.1:5;
        y1 = repmat(y(end),1,length(x1));

        x = [x, x1];
        y = [y, y1];

        t = pi/2:0.1:3*pi/2;
        x1 = 4.5.*cos(t) + 5;
        y1 = 4.5.*sin(t) + 5;

        x = [x, x1];
        y = [y, y1];

        d = -max(S.system.Lh(:,2))*S.config.Nt*log(tan(0.1/2));
        y1 = y(1)-d:0.1:y(1);
        x1 = x(1).*ones(size(y1));

        y_aux = [y1,y];
        S.path.coordinates = [x1,x;y_aux];  
        
    elseif strcmp(type,'rectangular_terrace') || strcmp(type,'complex1_terrace')
        if strcmp(type,'rectangular_terrace')
%            bag = rosbag('/home/kernighan/2022-07-19-15-32-59.bag');             
            x = repmat(8.2,1,101);
            y = 20:0.1:30;
            x = [x, 8.2:-0.1:1.0];
            y = [y, repmat(30, 1, ceil((8.2-1.0)/0.1)+1)];
            x = [x, repmat(0.8,1,101)];
            y = [y, 30:-0.1:20];
            x = [x, 0.8:0.1:8.2];
            y = [y, repmat(20, 1, ceil((8.2-0.8)/0.1)+1)];
           
            d = -max(S.system.Lh(:,2))*S.config.Nt*log(tan(0.1/2));
            y1 = y(1)-d:0.1:y(1);
            x1 = x(1).*ones(size(y1));

            y_aux = [y1,y];
            y_aux = y_aux+1.25;
            S.path.coordinates = [x1,x;y_aux];                      
           
        elseif strcmp(type,'complex1_terrace')
           bag = rosbag('/home/nahuel/Dropbox/PosDoc/AC3E/NMHE for incremental encoders/2022-07-19-15-36-09.bag'); 
%            bag = rosbag('/home/kernighan/Documents/mhe-mpc-for-N-trailers/Simulaciones/ACADO/2022-08-10-17-00-40.bag');
%            bag = rosbag('/home/nahuel/Dropbox/PosDoc/AC3E/NMHE-NMPC-for-N-trailer/mhe-mpc-for-N-trailers/Simulaciones/ACADO/2022-07-19-15-36-09.bag');            
%             bag         = rosbag('/home/nahuel/Dropbox/PosDoc/AC3E/NMHE for incremental encoders//test_betas_3.bag');            
            bSel        = select(bag, 'Topic', 'fix');
            msgStructs  = readMessages(bSel);

            lat         = NaN(1,bSel.NumMessages);
            lon         = NaN(1,bSel.NumMessages);
            xraw        = NaN(1,bSel.NumMessages);
            yraw        = NaN(1,bSel.NumMessages);
            xcor        = NaN(1,bSel.NumMessages);
            ycor        = NaN(1,bSel.NumMessages);
            % ESQUINA DE LA AZOTEA TOMDA COMO ORIGEN DE COORDENADAS
            lat0        = -33.03422;%S.ROS.LAT0; %   = -33.03422;
            lon0        = -71.591885;%S.ROS.LON0; %  = -71.591885;
           
            %
            sdpvar a11 a12 a21 a22;
            sdpvar x1bar y1bar x2bar y2bar;
            % Besides reference point, two more are needed to obtaint the local
            % reference frame
            % Coordinates of point (x, 0)
            lat_local_coord_1 = -33.03421;%S.ROS.local_coord_1.lat;
            long_local_coord_1 = -71.591825;%S.ROS.local_coord_1.long;

            [x1, y1]    = latlon2xy(lat_local_coord_1, long_local_coord_1, lat0, lon0);
            x1          = x1*1000;
            y1          = y1*1000;
            x           = norm([x1 y1]);
            % Coordinates of point (0, y)
            lat_local_coord_2 = -33.034158333;%S.ROS.local_coord_2.lat;
            long_local_coord_2 = -71.591905;%S.ROS.local_coord_2.long;

            [x2, y2]    = latlon2xy(lat_local_coord_2, long_local_coord_2, lat0, lon0);
            x2          = x2*1000;
            y2          = y2*1000;
            y           = norm([x2 y2]);
            %
            A           = [a11 a12; a21 a22];
            v           = [x1; y1; x2; y2];
            b           = [x1bar; y1bar; x2bar; y2bar];
            Constraints = [[A zeros(2); zeros(2) A]*v - b == zeros(4,1); x1bar*x2bar + y1bar*y2bar == 0];
            %
            obj         = (x1bar - x)^2 + (y2bar - y)^2;
            % 
            optimize(Constraints, obj);
            %
            Mtx   = value(A);

            %

            for i=1:bSel.NumMessages
                lat(i)  = msgStructs{i}.Latitude;
                lon(i)  = msgStructs{i}.Longitude;
                [xm,ym] = latlon2xy(lat(i),lon(i),lat0,lon0);
                xraw(i) = xm*1000;
                yraw(i) = ym*1000;
                xy      = Mtx*[xraw(i);yraw(i)];
                xcor(i) = xy(1);
                ycor(i) = xy(2);
            end

            x = smooth(xcor,35)';
            y = smooth(ycor,35)';

            S.path.coordinates = [x; y];            
        end
%         for i=1:15
%             xcor = smooth(xcor,'lowess')';
%             ycor = smooth(ycor,'lowess')';
%         end
        
    elseif strcmp(type,'circular_terrace')           
        a = 4;
        b = 4;
%         b = 0.4;
        t = 0:0.0005:2*pi;        
        x = a*cos(t);
        y = b*sin(t);
        %
        S.path.coorection_x = 5.5;
        S.path.coorection_y = 4.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = a*cos(S.path.s)+S.path.coorection_x;
        S.path.fy           = b*sin(S.path.s)+S.path.coorection_y;
        %
        S.path.st           = casadi.MX.sym('st');
        S.path.ftheta       = pi/2 + (1/S.config.Ts)*S.path.st * atan(S.path.ds);
        S.path.r            = [ S.path.fx;  
                                S.path.fy;
                                S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type,'ellipsoidal1_terrace')           
           y = 18:0.1:23.5;
           x = repmat(8,1,length(y));
           
           t = 0:0.01:0.75;
           x = [x, 4.5 + 3.5.*cos(2*pi*t)];
           y = [y, 23.5 + 2.5.*sin(2*pi*t)];
           
           x1 = x(end):0.05:5.5;
           y1 = repmat(y(end),1,length(x1));
           
           x = [x,x1];
           y = [y,y1];
           
           y1 = y(end):-0.05:y(1)-2;
           x1 = repmat(x(end),1,length(y1));
           
           x = [x,x1];
           y = [y,y1];
           
           t = pi:0.001:2*pi;
           
           x1 = 1.25*cos(t) + 6.75;
           y1 = 1.25*sin(t) + 16;
           
           x = [x,x1];
           y = [y,y1];
           
           y1 = y(end):0.05:y(1);
           x1 = repmat(x(end),1,length(y1));
           
           x = [x, x1];
           y = [y, y1];
           
           x = x-1;
           y = y+4;

           S.path.coordinates = [x; y];
    elseif strcmp(type,'ellipsoidal2_terrace')           
           y = 29:-0.1:23.5;
           x = repmat(8,1,length(y));
           
           t = 0:0.01:1;
           x = [x, 4.5 + 3.5.*cos(-2*pi*t)];
           y = [y, 23.5 + 2.5.*sin(-2*pi*t)];
           
           x = x;
           y = y+2;

           S.path.coordinates = [x; y];
    elseif strcmp(type,'agricultural')
            x0 = 5;
            y0 = 1;
            
            factor = 3;

            yinit = 0;
            
            d1 = 1;
            d2 = 0.5;
            d3 = 1;
            
            a1 = 2.5/factor;
            b1 = 2.5/factor;
            a2 = 3.5/factor;
            b2 = 2.5/factor;
            a3 = 2.5/factor;
            b4 = 2.5/factor;
            
            h  = 10/factor;
            
            delta = 0.00001;
            alpha1 = 25*pi/180;
            alpha2 = 35*pi/180;
            alpha3 = 25*pi/180;
            
            xc1 = x0+d1/2;
            yc1 = y0+h+b1;
            
            xc2 = x0+d2/2;
            yc2 = y0+b2;
            
            xc3 = x0+d1+d2+d3/2;
            yc3 = y0+h+b1;
            
            t1  = 3*pi/2-alpha1:-delta:-pi/2+alpha1;
            xe1 = xc1+a1*cos(t1);
            ye1 = yc1+b1*sin(t1);
            
            t2  = pi-(pi/2-alpha2):delta:2*pi+(pi/2-alpha2);
            xe2 = xc2+a2*cos(t2);
            xe2 = xe2 + xe1(end) - xe2(1);
            ye2 = yc2+b2*sin(t2);
            
            t3  = 3*pi/2-alpha3:-delta:-pi/2+alpha3;
            xe3 = xc3+a3*cos(t3);
            xe3 = xe3+xe2(end)-xe3(1); 
            ye3 = yc3+b4*sin(t3);
            
            y   = 2+h/2:delta:ye1(1);
            x   = repmat(xe1(1),1,length(y));
            
            path = [x, xe1;y, ye1];
            
            y = ye1(end):-delta:ye2(1);
            x = repmat(xe1(end),1,length(y));
            
            path = [path, [x,xe2;y, ye2]];
            
            y = ye2(end):delta:ye3(1);
            x = repmat(xe2(end),1,length(y));
            
            path = [path, [x,xe3;y, ye3]];
            
            y = ye3(end):-delta:yinit-delta;
            x = repmat(xe3(end),1,length(y));
            
            path = [path,[x;y]];
            
            x = path(1,end):-delta:path(1,1);
            y = repmat(path(2,end),1,length(x));
            
            path = [path,[x;y]];
            
            y = path(2,end):delta:path(2,1);
            x = repmat(path(1,end),1,length(y));
            
            path = [path,[x;y]];
            S.path.coordinates = path;
    else    
        % REAL TRAJECTORY #####################################################
        % First segment
        % Path ________________________________________________________________
        S.path.x0               = 5;
        S.path.y0               = 5;    
        S.path.R1               = 8;
        S.path.R2               = 6;
        S.path.R3               = 4;
        S.path.R4               = 3;
        S.path.R5               = 1;
        S.path.L1               = 20;
        S.path.L2               = 10;
        S.path.L3               = 2*(S.path.R1+S.path.R2+S.path.R3+S.path.R4+S.path.R5);
        S.path.deltaL           = 0.1;
        %
        S.path.coordinates = [repmat(S.path.x0, 1, S.path.L1/S.path.deltaL); S.path.y0:S.path.deltaL:S.path.L1+S.path.y0-S.path.deltaL];
        % First circunference
        N   = pi*S.path.R1/S.path.deltaL - S.path.deltaL;
        t   = pi:-1/N:0;
        xc1 = S.path.coordinates(1,end);
        yc1 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc1 + S.path.R1 + S.path.R1.*cos(t); yc1 + S.path.R1.*sin(t)]];
        % Second segment
        x2 = S.path.coordinates(1,end);
        y2 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x2, 1, S.path.L2/S.path.deltaL); y2:-S.path.deltaL:y2-S.path.L2+S.path.deltaL] ];
        % Second circunference
        N   = pi*S.path.R2/S.path.deltaL - S.path.deltaL;
        t   = pi:1/N:2*pi;
        xc2 = S.path.coordinates(1,end);
        yc2 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc2 + S.path.R2 + S.path.R2.*cos(t); yc2 + S.path.R2.*sin(t)]];
        % Third segment
        x3 = S.path.coordinates(1,end);
        y3 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x3, 1, S.path.L2/S.path.deltaL); y3:S.path.deltaL:y3+S.path.L2-S.path.deltaL] ];
        % Third circunference
        N   = pi*S.path.R3/S.path.deltaL - S.path.deltaL;
        t   = pi:-1/N:0;
        xc3 = S.path.coordinates(1,end);
        yc3 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc3 + S.path.R3 + S.path.R3.*cos(t); yc3 + S.path.R3.*sin(t)]];
        % Fourth segment
        x4 = S.path.coordinates(1,end);
        y4 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x4, 1, S.path.L2/S.path.deltaL); y4:-S.path.deltaL:y4-S.path.L2+S.path.deltaL] ];
        % Fourth circunference
        N   = pi*S.path.R4/S.path.deltaL - S.path.deltaL;
        t   = pi:1/N:2*pi;
        xc4 = S.path.coordinates(1,end);
        yc4 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc4 + S.path.R4 + S.path.R4.*cos(t); yc4 + S.path.R4.*sin(t)]];
        % Fifth segment
        x5 = S.path.coordinates(1,end);
        y5 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x5, 1, S.path.L2/S.path.deltaL); y5:S.path.deltaL:y5+S.path.L2-S.path.deltaL] ];
        % Fifth circunference
        N   = pi*S.path.R5/S.path.deltaL - S.path.deltaL;
        t   = pi:-1/N:0;
        xc5 = S.path.coordinates(1,end);
        yc5 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc5 + S.path.R5 + S.path.R5.*cos(t); yc5 + S.path.R5.*sin(t)]];
        % Sexth segment
        x6 = S.path.coordinates(1,end);
        y6 = S.path.coordinates(2,end);

        y6_aux = y6:-S.path.deltaL:S.path.y0-0.5;%+S.path.deltaL;
        %
        S.path.coordinates = [S.path.coordinates, [repmat(x6, 1, length(y6_aux)); y6_aux] ];
        % Merge with first point
        x8 = S.path.coordinates(1,end);
        y8 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [x8:-S.path.deltaL:S.path.x0; [repmat(y8, 1, (S.path.L3/S.path.deltaL)/2), repmat(S.path.y0, 1, (S.path.L3/S.path.deltaL)/2)]] ];
        %
        yval = S.path.coordinates(2,end);
        for i=length(S.path.coordinates):-1:2
            if S.path.coordinates(2,i) ~= yval
                S.path.coordinates(1,i+1) = S.path.coordinates(1,i);
                break;
            end
        end
        S.path.coordinates = [S.path.coordinates, S.path.coordinates(:,1)];
        %
        S.path.num_points_ttaj = length(S.path.coordinates);
%         % Gen path to intial point x0
%         x9 = S.path.coordinates(1,end);
%         y9 = S.path.coordinates(2,end);
%         trajX = x9-S.path.deltaL:-S.path.deltaL:S.init_condition.x0(2*S.config.Nt+2)-2*(1+S.config.Nt)*(S.system.Lh1+S.system.L1);
%         trajY = y9.*ones(size(trajX));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
%         %
%         x10 = S.path.coordinates(1,end);
%         y10 = S.path.coordinates(2,end);    
%         trajY = y10-S.path.deltaL:-S.path.deltaL:S.init_condition.x0(2*S.config.Nt+3);
%         trajX = x10.*ones(size(trajY));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
%         %
%         x11 = S.path.coordinates(1,end);
%         y12 = S.path.coordinates(2,end);
%         trajX = x11+S.path.deltaL:S.path.deltaL:S.init_condition.x0(2*S.config.Nt+2);
%         trajY = y12.*ones(size(trajX));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
        % Remove inconsistencies
        pos_to_del = [];
        for i=2:length(S.path.coordinates)
            dx = S.path.coordinates(1,i)-S.path.coordinates(1,i-1);
            dy = S.path.coordinates(2,i)-S.path.coordinates(2,i-1);
            if dx==0 && dy ==0
                pos_to_del = [pos_to_del,i];                
            end
        end
        S.path.coordinates(:,pos_to_del) = [];
    end
    % Compute length of the path ------------------------------------------
    len = 0;
    for i=2:length(S.path.coordinates)
        len = len + norm(S.path.coordinates(:,i)-S.path.coordinates(:,i-1));
    end
    S.path.length = len;
    % Gen Slippage vector for each coordinate of the path -----------------
%     S.path.slippage = ones(S.system.np,length(S.path.coordinates));
end

function S = genObstacle(S,x,y,r)
    % It genearte on osbtacle. Other functions are inc hange of detecting
    % which ones of them are static and dynamics
    S.path.listOfObs = [S.path.listOfObs; [x, y, r]];
end

function S = gen_dynamic_and_solvers(S)
    % System's dimension __________________________________________________        
    S.system.nu      = 2;                        % input dimension: [w0,v0]
    S.system.ny      = S.config.Nt + 3;
    S.system.nv      = S.system.ny;  
    S.system.nq      = 2*S.config.Nt + 1 + 2*(S.config.Nt+1);
    % System's parameters -------------------------------------------------
    S.system.Lh1     = 0.342;
    S.system.L1      = 1.08;
    S.system.Lh2     = 0;
    S.system.L2      = 0.78;
    % ---------------------------------------------------------------------
    S.system.Lh3     = 0.15;%-0.342;
    S.system.L3      = 0.7;
    S.system.Lh4     = -0.15;%0.342;
    S.system.L4      = 0.7;
    S.system.Lh5     = 0;
    S.system.L5      = 0.7;
    S.system.Lh6     = -0.342;
    S.system.L6      = 0.7;
    S.system.Lh7     = 0.342;
    S.system.L7      = 0.7;
    S.system.Lh8     = 0;
    S.system.L8      = 0.7;
    S.system.Lh9     = 0;%0.3;%0.342;%0.048;
    S.system.L9      = 0.25;%0.78;%0.229;
    S.system.Lh10    = 0;%0.3;%0.048;
    S.system.L10     = 0.25;%0.78;%0.229;
    %
    S.system.Lhi     = [S.system.Lh1;S.system.Lh2;S.system.Lh3;S.system.Lh4;S.system.Lh5;S.system.Lh6;S.system.Lh7;S.system.Lh8;S.system.Lh9;S.system.Lh10]; 
    S.system.Li      = [S.system.L1;S.system.L2;S.system.L3;S.system.L4;S.system.L5;S.system.L6;S.system.L7;S.system.L8;S.system.L9;S.system.L10]; 
    %
    S.system.r       = 0.15;           % wheel radius
    S.system.width   = 0.67;%0.229;    % vehicles's width    
    S.system.long    = 0.99;%0.229;    % vehicles's long
    % For plotting purposes
    S.system.XYtracAxe           = [-S.system.width/2 S.system.width/2; 0 0];
    S.system.XYtracWheelLeft     = [-S.system.width/2 -S.system.width/2; S.system.r/2 -S.system.r/2];
    S.system.XYtracWheelRight    = [S.system.width/2 S.system.width/2; -S.system.r/2 S.system.r/2];
%     S.system.XYtracBody          = [-S.system.width/2*0.75 0 S.system.width/2 -S.system.width/2*0.75; -S.system.width/2*0.75 S.system.width/2*2.75 -S.system.width/2*0.75 -S.system.width/2*0.75];

    VerticesHusky                = [-S.system.width/2  -S.system.long/2; -S.system.width/2 S.system.long/2; S.system.width/2 S.system.long/2; S.system.width/2  -S.system.long/2];
    S.system.XYtracBody          = Polyhedron(VerticesHusky);

    %
    S.system.XYtrailerAxe        = [-S.system.width/2 S.system.width/2; 0 0];
    S.system.XYtrailerWheelLeft  = [-S.system.width/2 -S.system.width/2; S.system.r/2 -S.system.r/2];
    S.system.XYtrailerWheelRight = [S.system.width/2 S.system.width/2; -S.system.r/2 S.system.r/2];
S.system.XYtrailerLoad       = {};
    S.system.XYtrailerLongAxe    = [];
    for i=1:S.config.Nt
        S.system.XYtrailerLongAxe = [S.system.XYtrailerLongAxe; [0 0; 0 S.system.Li(i)]];

%         [S.system.XYtrailerLoad; [-S.system.width/2*0.9 -S.system.width/2*0.9 S.system.width/2*0.9 S.system.width/2*0.9 -S.system.width/2*0.9; -S.system.r/2*1.2 S.system.Li(i)*0.65 S.system.Li(i)*0.65 -S.system.r/2*1.2 -S.system.r/2*1.2]];
S.system.XYtrailerLoad{i} = Polyhedron([-S.system.width/2 -S.system.r/2; -S.system.width/2 (S.system.Li(i)-S.system.r/2)*0.6; S.system.width/2 (S.system.Li(i)-S.system.r/2)*0.6; S.system.width/2 -S.system.r/2]);
    end
    % Casadi variabes ----------------------------------------------------
    S.dynamic.q      = casadi.MX.sym('q',S.system.nq);
    S.dynamic.u      = casadi.MX.sym('u',S.system.nu);
    % System's variables  -------------------------------------------------      
    switch(S.config.Nt)
        case 1            
            % CASADI FORMULATION ******************************************
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1, sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];
            S.system.J       = [S.system.J1];
            %
            S.dynamic.f_rhs  = [ [1 0] * (eye(2)-S.system.J1) * S.dynamic.u;...
                                 %
                                 [1 0] * S.dynamic.u;...
                                 [1 0] * S.system.J1 * S.dynamic.u;...
                                 %
                                 [0 cos(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                 [0 sin(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                 [0 cos(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u;...
                                 [0 sin(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u  ];
            %
            S.system.Lh      = [S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 2
            % CASADI FORMULATION ******************************************
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2, sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)), cos(S.dynamic.q(2))];            
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1, sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];            
            S.system.J       = [S.system.J2, S.system.J1];

            S.dynamic.f_rhs  = [    [1 0] * (eye(2)-S.system.J1) * S.dynamic.u;
                                    [1 0] * (eye(2)-S.system.J2) * S.system.J1 * S.dynamic.u;
                                    %
                                    [1 0] * S.dynamic.u;
                                    [1 0] * S.system.J1 * S.dynamic.u;
                                    [1 0] * S.system.J2 * S.system.J1 * S.dynamic.u;
                                    %
                                    [0 cos(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.Nt+3))] * S.system.J2 * S.system.J1 * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.Nt+3))] * S.system.J2 * S.system.J1 * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 3
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+3))] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))] * S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.Nt+1))] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.Nt+1))] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 4
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))] * S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+3))] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))] * S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+4))] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+4))] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.Nt+1))] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.Nt+1))] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh4,S.system.L4;...
                                    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 5
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))] *S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))] *S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+3))] *S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))] *S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+4))] *S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+4))] *S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+5))] *S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+5))] *S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.Nt+1))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.Nt+1))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh5,S.system.L5;...
                                    S.system.Lh4,S.system.L4;...
                                    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 6
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))]*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))]*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh6,S.system.L6;...
                                    S.system.Lh5,S.system.L5;...
                                    S.system.Lh4,S.system.L4;...
                                    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 7
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1)*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh7,S.system.L7;...
                                    S.system.Lh6,S.system.L6;...
                                    S.system.Lh5,S.system.L5;...
                                    S.system.Lh4,S.system.L4;...
                                    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
        case 8
            S.system.J8      = [-S.system.Lh8 * cos(S.dynamic.q(8)) / S.system.L8 sin(S.dynamic.q(8))/S.system.L8; S.system.Lh8*sin(S.dynamic.q(8)) cos(S.dynamic.q(8))];
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J8, S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1)*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0]*S.system.J1* S.dynamic.u;...
                                   [1 0]*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+9))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+9))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh8,S.system.L8;...
                                    S.system.Lh7,S.system.L7;...
                                    S.system.Lh6,S.system.L6;...
                                    S.system.Lh5,S.system.L5;...
                                    S.system.Lh4,S.system.L4;...
                                    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
        case 9
            S.system.J9      = [-S.system.Lh9 * cos(S.dynamic.q(9)) / S.system.L9 sin(S.dynamic.q(9))/S.system.L9; S.system.Lh9*sin(S.dynamic.q(9)) cos(S.dynamic.q(9))];
            S.system.J8      = [-S.system.Lh8 * cos(S.dynamic.q(8)) / S.system.L8 sin(S.dynamic.q(8))/S.system.L8; S.system.Lh8*sin(S.dynamic.q(8)) cos(S.dynamic.q(8))];
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J9, S.system.J8, S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1)*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J9)*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+9))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+9))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+10))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+10))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh9,S.system.L9;...
                                    S.system.Lh8,S.system.L8;...
                                    S.system.Lh7,S.system.L7;...
                                    S.system.Lh6,S.system.L6;...
                                    S.system.Lh5,S.system.L5;...
                                    S.system.Lh4,S.system.L4;...
                                    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
        case 10
            S.system.J10     = [-S.system.Lh10 * cos(S.dynamic.q(10)) / S.system.L10 sin(S.dynamic.q(10))/S.system.L10; S.system.Lh10*sin(S.dynamic.q(10)) cos(S.dynamic.q(10))];
            S.system.J9      = [-S.system.Lh9 * cos(S.dynamic.q(9)) / S.system.L9 sin(S.dynamic.q(9))/S.system.L9; S.system.Lh9*sin(S.dynamic.q(9)) cos(S.dynamic.q(9))];
            S.system.J8      = [-S.system.Lh8 * cos(S.dynamic.q(8)) / S.system.L8 sin(S.dynamic.q(8))/S.system.L8; S.system.Lh8*sin(S.dynamic.q(8)) cos(S.dynamic.q(8))];
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J10, S.system.J9, S.system.J8, S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];            
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1)*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J9)*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J10)*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+1))]*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+2))]*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+3))]*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+4))]*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+9))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+9))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+10))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+10))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.Nt+11))]*S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.Nt+11))]*S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1*S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh10,S.system.L10;...
                                    S.system.Lh9,S.system.L9;...
                                    S.system.Lh8,S.system.L8;...
                                    S.system.Lh7,S.system.L7;...
                                    S.system.Lh6,S.system.L6;...
                                    S.system.Lh5,S.system.L5;...
                                    S.system.Lh4,S.system.L4;...
                                    S.system.Lh3,S.system.L3;...
                                    S.system.Lh2,S.system.L2;...
                                    S.system.Lh1,S.system.L1];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
    end
    % Generate CASADI functions and integrators ***************************
    S.dynamic.f     = casadi.Function('f_rhs', {S.dynamic.q,S.dynamic.u}, {S.dynamic.f_rhs});    
    opts            = struct('main',true,'mex',true);
    S.dynamic.f.generate('f.c',opts)
    mex f.c -largeArrayDims;
    % Integrator ----------------------------------------------------------
    if strcmp(S.config.integrator,'RK4')
        k1              = S.dynamic.f(S.dynamic.q, S.dynamic.u);
        k2              = S.dynamic.f(S.dynamic.q + S.config.Ts / 2 * k1, S.dynamic.u);
        k3              = S.dynamic.f(S.dynamic.q + S.config.Ts / 2 * k2, S.dynamic.u);
        k4              = S.dynamic.f(S.dynamic.q + S.config.Ts * k3, S.dynamic.u);
        x_rk4           = S.dynamic.q + S.config.Ts / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    elseif strcmp(S.config.integrator,'euler')
        x_rk4           = S.dynamic.q + S.config.Ts*S.dynamic.f(S.dynamic.q, S.dynamic.u);
    end
    S.dynamic.FNt   = casadi.Function('FNt', {S.dynamic.q, S.dynamic.u}, {x_rk4});
    S.dynamic.FNt.generate('FNt.c',opts);
    mex FNt.c -largeArrayDims;
    % Output of the system ------------------------------------------------
    S.dynamic.h_rhs  = S.dynamic.q(S.config.outputs);
    S.dynamic.h      = casadi.Function('h', {S.dynamic.q}, {S.dynamic.h_rhs});
end

function S = init_mheCasadi(S)
    S.mpc.box_constraints               = struct;
    S.mpc.box_constraints.QluBounds     = [];
    S.mpc.box_constraints.WluBounds     = [];
    S.mpc.box_constraints.VluBounds     = [];
    S.mpc.box_constraints.ZluBounds     = [repmat(-0.01,S.config.Nt,1), repmat(0.01,S.config.Nt,1)];
    S.mpc.box_constraints.UluBounds     = [];
    %
    S.mpc.Mtxs                          = struct;
    
    if S.config.SIM == false % G2T
        S.mpc.Mtxs.Q    = diag([0.5.*ones(1,S.config.Nt),1,0.25.*ones(1,S.config.Nt),ones(1,2),0.25.*ones(1,2*S.config.Nt)]);%,0.2,0.2]);
        S.mpc.Mtxs.R    = diag([0.25.*ones(1,S.config.Nt),0.65,1,1]);%,1,1]);
        S.mpc.Mtxs.P    = 1e6.*eye(S.system.nq);
        S.mpc.Mtxs.Z    = 1e0.*eye(S.config.Nt);
        S.mpc.Mtxs.dU   = eye(S.system.nu);
    else
        S.mpc.Mtxs.Q    = diag([0.5.*ones(1,S.config.Nt),1,0.25.*ones(1,S.config.Nt),ones(1,2),0.25.*ones(1,2*S.config.Nt)]);%,0.2,0.2]);
        S.mpc.Mtxs.R    = diag([0.25.*ones(1,S.config.Nt),0.65,1,1]);%,1,1]);
        S.mpc.Mtxs.P    = 1e6.*eye(S.system.nq);
        S.mpc.Mtxs.Z    = 1e0.*eye(S.config.Nt);
        S.mpc.Mtxs.dU   = eye(S.system.nu);
    end
    % nlmheCasadiNt(Ne,x,u,Nt,f_rhs,h_rhs,Mtxs,nq,nu,ny,boxConst,q0bar,Ts,dimensions)
    S.algorithms.mheCasadi = nlmheCasadiNt(S.config.Ne,S.dynamic.q,S.dynamic.u,S.config.Nt,S.dynamic.f_rhs,S.dynamic.h_rhs,...
        S.mpc.Mtxs,S.system.nq,S.system.nu,S.system.ny,S.mpc.box_constraints,S.init_condition.x0bar,S.config.Ts,S.system.Lh);
    %
    setSigma(S.algorithms.mheCasadi, 1);
    setC(S.algorithms.mheCasadi, 1e7);
end

function S = fill_mhe(S)
    for imhe=1:S.config.Ne
        updateMeasurement(S.algorithms.mheCasadi, full(S.dynamic.h(S.init_condition.x0bar)));
        updateInput(S.algorithms.mheCasadi,zeros(S.system.nu,1));       
    end
    updateMeasurement(S.algorithms.mheCasadi, full(S.dynamic.h(S.init_condition.x0bar)));
end

function S = init_fnmppc(S)
    S.fnmppc.optiMPC = casadi.Opti();
    S.fnmppc.qMPC    = S.fnmppc.optiMPC.variable(S.system.nq,S.config.Nc+1);
    S.fnmppc.u       = S.fnmppc.optiMPC.variable(S.system.nu,S.config.Nc);
    S.fnmppc.qref    = S.fnmppc.optiMPC.parameter(3);
    S.fnmppc.qinit   = S.fnmppc.optiMPC.parameter(S.system.nq);
    S.fnmppc.u_last  = S.fnmppc.optiMPC.parameter(S.system.nu);
    S.fnmppc.ObsSD    = S.fnmppc.optiMPC.parameter(S.config.totNumObs,5); % [x,y,r,delta x, delta y]' x Nd: with delta x and delta y, I can predict the future position by a linear extrapolation

    S.fnmppc.Qmpc    = diag([1 10 10]);
    S.fnmppc.QNcmpc  = S.fnmppc.Qmpc;
    S.fnmppc.Rmpc    = diag([0.05 0.1]);
    S.fnmppc.Raccel  = diag([1 1]);

    S.fnmppc.UluBounds  = [-2 2;...
                           -1 1];
    S.fnmppc.dUluBounds = S.config.Ts.*[-6 6;...  % (rad/s^2), data from: https://answers.ros.org/question/335763/what-is-the-max-linear-acceleration-of-husky-a200/
                                        -3 3];    % (m/s^2)
    S.fnmppc.QluBounds  = [repmat([-110*pi/180 110*pi/180],S.config.Nt,1); repmat([-inf inf],S.config.Nt+1,1); repmat([-inf inf],2*(S.config.Nt+1),1)];

    Jmpc = 0;
    for k=1:S.config.Nc+1        
        if k < S.config.Nc+1
            Jmpc = Jmpc + (S.fnmppc.qMPC([S.config.Nt+1,2*S.config.Nt+2:2*S.config.Nt+3],k)-S.fnmppc.qref).'* S.fnmppc.Qmpc *(S.fnmppc.qMPC([S.config.Nt+1,2*S.config.Nt+2:2*S.config.Nt+3],k)-S.fnmppc.qref);
            %
            S.fnmppc.optiMPC.subject_to( S.fnmppc.qMPC(:,k+1) == S.dynamic.FNt(S.fnmppc.qMPC(:,k),S.fnmppc.u(:,k)) );           
            if k==1
                Jmpc = Jmpc + (S.fnmppc.u(:,k)-S.fnmppc.u_last)'*S.fnmppc.Raccel*(S.fnmppc.u(:,k)-S.fnmppc.u_last);
            else
                Jmpc = Jmpc + (S.fnmppc.u(:,k)-S.fnmppc.u(:,k-1))'*S.fnmppc.Raccel*(S.fnmppc.u(:,k)-S.fnmppc.u(:,k-1));
            end        
            %
            Jmpc = Jmpc + S.fnmppc.u(:,k).'* S.fnmppc.Rmpc *S.fnmppc.u(:,k);
        else
            Jmpc = Jmpc + (S.fnmppc.qMPC([S.config.Nt+1,2*S.config.Nt+2:2*S.config.Nt+3],k)-S.fnmppc.qref).'* S.fnmppc.QNcmpc *(S.fnmppc.qMPC([S.config.Nt+1,2*S.config.Nt+2:2*S.config.Nt+3],k)-S.fnmppc.qref);
        end
        % handle dynamis obstacles
        if k>1
            for j=1:S.config.totNumObs
                for i=2:S.config.Nt+1
                    for l=0:1/4:1
                        vx = S.fnmppc.qMPC(2*S.config.Nt+1+(i-1)*2+1,k) - S.fnmppc.qMPC(2*S.config.Nt+1+(i-2)*2+1,k);
                        vy = S.fnmppc.qMPC(2*S.config.Nt+1+i*2,k) - S.fnmppc.qMPC(2*S.config.Nt+1+(i-1)*2,k);
   
                        S.fnmppc.optiMPC.subject_to( S.fnmppc.ObsSD(j,3)^2 < (S.fnmppc.qMPC(2*S.config.Nt+1+(i-2)*2+1,k)+vx*l - (S.fnmppc.ObsSD(j,1)+S.fnmppc.ObsSD(j,4)*(k-1)*S.config.Ts))^2 +... 
                                                                                 (S.fnmppc.qMPC(2*S.config.Nt+1+(i-1)*2,k)+vy*l - (S.fnmppc.ObsSD(j,2)+S.fnmppc.ObsSD(j,5)*(k-1)*S.config.Ts))^2);
                    end
                end
            end
        end
    end
    %
    S.fnmppc.optiMPC.minimize(Jmpc);
    %
    S.fnmppc.optiMPC.subject_to(S.fnmppc.UluBounds(:,1) < S.fnmppc.u < S.fnmppc.UluBounds(:,2));
    S.fnmppc.optiMPC.subject_to( S.fnmppc.dUluBounds(:,1) <= S.fnmppc.u(:,1)-S.fnmppc.u_last <= S.fnmppc.dUluBounds(:,2) );
    S.fnmppc.optiMPC.subject_to( S.fnmppc.dUluBounds(:,1) <= S.fnmppc.u(:,2:end)-S.fnmppc.u(:,1:end-1) <= S.fnmppc.dUluBounds(:,2) );
    S.fnmppc.optiMPC.subject_to(S.fnmppc.qMPC(:,1) == S.fnmppc.qinit);
    S.fnmppc.optiMPC.subject_to(S.fnmppc.QluBounds(:,1) < S.fnmppc.qMPC(:,2:end) < S.fnmppc.QluBounds(:,2));
    %
    p_optsMPC = struct('expand',true);
    s_optsMPC = struct('sb','yes','print_level',0,'gamma_theta',1e-1,'jacobian_approximation','exact','fast_step_computation','yes','warm_start_init_point','yes'); % 
    % Here effectively define the solver: every time this method is called, the
    % MPC probelm is solved
    S.fnmppc.optiMPC.solver('ipopt',p_optsMPC,s_optsMPC);
    S.fnmppc.optiMPC.set_value(S.fnmppc.u_last,[0;0]);

    S.mpc.Controls = zeros(S.system.nu,1); % Use same vector for all algorithms
    % Init with no obsatcles
    S.fnmppc.optiMPC.set_value(S.fnmppc.ObsSD,zeros(S.config.totNumObs,5)); 
end

function S = init_mpc(S)
    % PARAMETERS MPC
    S.mpc.box_constraints               = struct;    
    S.mpc.box_constraints.QluBounds     = [repmat([-120*pi/180 120*pi/180],S.config.Nt,1); repmat([-inf inf],S.config.Nt+1,1); repmat([-inf inf],2*(S.config.Nt+1),1)];%; [-1 1; -5 5]];
    if ~S.config.stab.compFlg
        S.mpc.box_constraints.QNluBounds    = [repmat([-120*pi/180 120*pi/180],S.config.Nt,1); repmat([-inf inf],S.config.Nt+1,1); repmat([-inf inf],2*(S.config.Nt+1),1)];%; [-1 1; -5 5]];
    else
        S.mpc.box_constraints.QNluBounds    = [repmat([-110*pi/180 110*pi/180],S.config.Nt,1); [S.config.stab.thetaMin S.config.stab.thetaMax] ;repmat([-inf inf],S.config.Nt,1); ...
                                              [S.config.stab.xMin S.config.stab.xMax]; [S.config.stab.yMin S.config.stab.yMax]; repmat([-inf inf],2*(S.config.Nt),1)];%; [-1 1; -5 5]];
    end


    S.mpc.box_constraints.UluBounds     = [-2 2;...
                                           -1 1];
    S.mpc.box_constraints.dUluBounds    = S.config.Ts.*[-6 6;...  % (rad/s^2), data from: https://answers.ros.org/question/335763/what-is-the-max-linear-acceleration-of-husky-a200/
                                                        -3 3];    % (m/s^2)
    %
    S.mpc.Mtxs                          = struct;   
    S.mpc.Mtxs.Q                        = diag([1 10 10]);
    S.mpc.Mtxs.QN                       = S.mpc.Mtxs.Q;
    S.mpc.Mtxs.R                        = diag([0.05 0.1]);
    %
    S.mpc.Mtxs.c                        = 0.001;          % loss function coefficient
    S.mpc.Mtxs.lossAmp                  = 10;            % loss function amplitude
    S.mpc.Mtxs.weightRect               = 3;            % penalisation of the distance to the line
    S.mpc.Mtxs.ampGauss                 = 100;%60;
    S.mpc.Mtxs.ampGaussDyn              = 100;%280;
    S.mpc.Mtxs.compContSet              = S.config.stab.compFlg;
    %
    S.mpc.Mtxs.Polyhedrons.Omega_O      = S.config.Omega_O;
    %
    S.mpc.dims                          = struct;
    S.mpc.dims.nq                       = S.system.nq;
    S.mpc.dims.nu                       = S.system.nu;
    %
    S.mpc.obsStrategy                   = 'gauss';
    %
%     S.mpc.mpcCasadi                     = nlmpcCasadiNt2(S.config.Nc,S.config.Nt,S.config.segmentTosteer,S.dynamic.FNt,S.mpc.Mtxs,S.mpc.dims,S.mpc.box_constraints,S.config.Ts,...
%                                             S.mpc.obsStrategy,S.config.numStObsSgmnt*(S.config.Nt+1),S.config.numMvObsSgmnt);

    S.mpc.mpcCasadi                     = nlmpcCasadiOptiNt2(S.config.Nc,S.config.Np,S.config.Nt,S.config.segmentTosteer,S.dynamic.FNt,S.mpc.Mtxs,S.mpc.dims,S.mpc.box_constraints,S.config.Ts,...
                                            S.mpc.obsStrategy,S.config.maxStaticObs*(S.config.Nt+1),S.config.maxMovingObs);
    %
    S.mpc.Controls                      = zeros(S.system.nu,1);
    %
    S.mpc.Controls_w0                   = 0;
    S.mpc.Controls_v0                   = 0;
    %
    S.mpc.Husky_w0                      = 0;
    S.mpc.Husky_v0                      = 0;
    % Initi stability parameters:
    Lmax                                = 7.847;
    Lmin                                = 0.707;
    rStab                               = 0.1;
    setLmax(S.mpc.mpcCasadi, Lmax);
    setLmin(S.mpc.mpcCasadi, Lmin);
    setrStab(S.mpc.mpcCasadi, rStab);
    computeCStab(S.mpc.mpcCasadi);
end

function S = init()
    % init fields of the structure
    S                 = struct;
    S.config          = struct;
    S.algorithms      = struct;
    S.noises          = struct;
    S.init_condition  = struct;
    S.path            = struct;
    S.obstacles       = struct;
    S.data            = struct;    
    S.plots           = struct;
    S.acado           = struct;
    S.acado.system    = struct;
    S.acado.states    = struct;
    S.acado.controls  = struct;
    S.exec_time       = struct;
    %
    S                 = reserve_notemp_memory(S);
    S                 = build_setup(S);
    S                 = gen_dynamic_and_solvers(S);
    S                 = makeSolverFindAltPath(S);
    S                 = makeSolverFindCoeffEllipse(S);
    S                 = isTargetOccludedMakeSolver(S);
    S                 = makeSolverFindReachablePoint(S);
    %
    S.path.path       ='infinity';% 'rectangular';
    S.path.name       = [S.path.path,'-N=',num2str(S.config.Nt),'-Ne=',num2str(S.config.Ne),'-Nc=',num2str(S.config.Nc),'-Ts=',num2str(S.config.Ts),'.mat'];
    S                 = gen_path(S,S.path.path);
    S.config.reference  = false;        
    len                 = length(S.path.coordinates);
    S.path.ref          = [zeros(S.config.Nt,len); zeros(S.config.Nt+1,len); S.path.coordinates; S.path.coordinates];

    S.path.listOfObs    = [];
    S.path.dynamicObs   = [];
    S.path.staticObs    = [];
    S.path.listOfObsStr = {};
    
    S.path.vehicleDims  = max(S.system.Lh(:,2))/2;
    S.path.safeMargin   = 0.3;%max(S.system.Lh(:,2))/2;
    S.path.radii        = [0.15 0.15 0.25 0.25 0.35 0.15 0.20 0.20 0.35 0.1 0.05 0.05 0.05]';    

    if strcmp(S.path.path,'infinity')
        S.path.xyObs    = [ 10.1,  8.1; 
                            8,   8;
                            6,   4;
                            4,   4;
                            2,   0.5;
                            0,   4;
                            1,   8;
                            4,   8;
                            5.2, 4.8;
                            7,   3;
                            8,   0;
                            9,   1;
                            11,  2.25]; % S.config.maxMovingObs+S.config.maxStaticObs must not exceed the total number of obstacles
        if S.config.SIM
%             S = genObstacle(S,1e6, 1e6, 0.35);
            if S.config.obsRand
                for i=1:S.config.totNumObs
                    S = genObstacle(S,rand*10,rand*10, 0.25+rand*0.5);
                end
            else
                for i=1:S.config.totNumObs
                    S = genObstacle(S,S.path.xyObs(i,1),S.path.xyObs(i,2),S.path.radii(i)+S.path.vehicleDims+S.path.safeMargin);
                end
            end
            S.path.listOfObsStr{1} = S.path.listOfObs;
        else
            S = genObstacle(S,1e6, 1e6, 0.35);
        end
    elseif strcmp(S.path.path,'flat_infinity')        
        if S.config.SIM
            S = genObstacle(S,1e6, 1e6, 0.45);
        else
            S = genObstacle(S,1e6, 1e6, 0.35);
        end
    elseif strcmp(S.path.path,'segment') || strcmp(S.path.path,'invert_segment')
        if S.config.SIM
            S = genObstacle(S,3, 6, 0.35);
            % S = genObstacle(S,5, 6, 0.5);
            S = genObstacle(S,8, 6, 0.45);
        else
            S = genObstacle(S,1e6, 1e6, 0.35);
        end        
    elseif strcmp(S.path.path,'rectangular')
        if S.config.SIM

        else           
            S = genObstacle(S,1e6, 1e6, 0.35);
        end
%         S = genObstacle(S,10,3, 0.25);
    elseif strcmp(S.path.path,'test')

    elseif strcmp(S.path.path,'monaco')

    elseif strcmp(S.path.path,'test_betas')

    elseif strcmp(S.path.path,'rectangular_terrace')        
        if ~S.config.reference
            S = genObstacle(S,1e6, 1e6, 0.35);
        else
            S = genObstacle(S,1e6, 1e6, 0.25);
%             S = genObstacle(S,6, 31.75, 0.25);
%             S = genObstacle(S,0.8,30.75, 0.25);
%             S = genObstacle(S,0.8,25.75, 0.25);
        end
    elseif strcmp(S.path.path,'complex1_terrace')
        
    elseif strcmp(S.path.path,'circular_terrace')           
        % S = genObstacle(S,1e6, 1e6, 0.35);
        S = genObstacle(S,5.5, 8.5, 0.35);
        S = genObstacle(S,1.5, 4.5, 0.35);
    elseif strcmp(S.path.path,'ellipsoidal1_terrace')           

    elseif strcmp(S.path.path,'ellipsoidal2_terrace')           

    else  
    end
end

function S = gen_init_conditions(S)
    S   = gen_x0(S); 
    S   = gen_x0bar(S);
%     S   = init_mheCasadi(S);
%     S   = init_ekfCasadi(S);
%     S   = fill_mhe(S);
    S   = init_mpc(S);
    S   = init_fnmppc(S);
%     S   = init_Michalek2017_robust_nSNT(S);
end

function S = init_flags_and_counters(S)
%
    S.path.reach_end_mhempc       = false;
    %
    S.config.time     = 0;
    S.config.iters    = 0;
end

function S = reserve_notemp_memory(S)
    S.data.mhempc.performance                       = struct;
    S.data.mhempc.performance.Psi_e      = [];
    S.data.mhempc.performance.Psi_e_vec  = [];
    S.data.mhempc.performance.maxPsi_e   = [];
    S.data.mhempc.performance.Psi_u      = [];    
    S.data.mhempc.performance.Psi_u_vec  = [];    
    S.data.mhempc.performance.maxPsi_u   = [];    
    S.data.mhempc.performance.time_tot   = [];
    S.data.mhempc.performance.est_err    = [];
    %
    S.data.num_sims                      = 1;
        %
    S.initConditionsStab                     = [];
    S.trajStab                               = {};
    S.UStab                                  = {};
    S.initConditionsStab                     = {};
    S.xNcStab                                = [];
end

function S = reserve_temp_memory(S)
    % Variables: states, outputs, etc. ____________________________________
    S.path.acumIndx               = 1;
    %
    S.data.xsim                   = [];
    S.data.ysim                   = [];
    S.sensors.velocities          = [];
    S.data.references             = [];
    %
    S.data.slip                   = [];
    S.data.procDist               = [];
    %
    S.data.measNoise              = [];
    S.data.UmeasNoise             = [];
    %
    S.data.test                   = [];
    %                        
    S.data.xest                   = [];
    %                        
    S.data.time_mhempc            = [];
    %
    S.path.references_mhempc      = [];%S.path.tgt;
    S.path.reach_end_mhempc       = false;
    S.path.last_tgt               = 0;
    %   
    S.exec_time.t_mhe             = [];
    S.exec_time.t_pfa             = [];
    S.exec_time.t_mpc             = [];
    S.exec_time.t_fnmppc          = [];
    S.exec_time.t_sensors         = [];
    S.exec_time.t_obsdetector     = [];
    S.exec_time.t_ctrl            = 0;
    S.exec_time.t_tot             = 0;
    S.exec_time.t_acum            = 0;
    S.exec_time.t_mis             = 0;
    S.exec_time.t_plt             = 0;
    %
    S.sensors.theta0              = [];
    S.sensors.theta2              = [];
    S.sensors.theta2              = [];
    S.sensors.vlp16               = {};
    % Sensor data from bolie phones
    S.mobile.data                 = {};
    %
    S.controller.ref.xyth           = [];
    S.controller.ref.a            = [];
    S.controller.ref.b            = [];
    S.controller.ref.c            = [];
    %
    S.controller.angle0_old       = 0;
    S.controller.modVelFac        = 1;
    %
    S.path.nearObs                = {};
    S.path.obsToMPC               = [];
%     S.path.dynamicObs             = [];
    S.path.dynObsAux              = [];
    S.path.dynObsXYcoords         = [];
    S.path.dynObsEllipsePrms.a               = [];
    S.path.dynObsEllipsePrms.b               = [];
    S.path.dynObsEllipsePrms.xc              = [];
    S.path.dynObsEllipsePrms.yc              = [];
    S.path.dynObsEllipsePrms.t0              = [];
    S.path.dynObsEllipsePrms.kt              = [];
    S.path.dynObsEllipsePrms.theta           = [];
    S.path.dynObsEllipsePrms.x0              = [];
    S.path.dynObsEllipsePrms.t               = [];
    S.path.dynObsEllipsePrms.tOld            = [];
    S.path.dynObsEllipsePrmsAux.x0           = [];
    S.path.dynObsEllipsePrmsAux.tOld         = 0;
    %
    S.path.posObs                            = {};
    %
    S.path.occluded.x0                       = [];
    S.path.occluded.flg                      = false;
    S.path.occluded.p1                       = [];
    S.path.occluded.p2                       = [];
    S.path.occluded.p3                       = [];
    %
    S.stability.q0Opt                        = [];
    S.controller.ref.xReachable              = [];
    S.controller.ref.yReachable              = [];
    S.controller.ref.thetaReachable          = [];
    %
    if S.config.video && ~S.config.stab.compFlg
        nameVid = strcat(['obsAvoidanceGNT-Nt',num2str(S.config.Nt),'-Nc',num2str(S.config.Nc),'-Np',num2str(S.config.Np),'-numMovObs',num2str(S.config.maxMovingObs),'-numStatObs',num2str(S.config.maxStaticObs),'-mthd-',S.config.method,'.avi']);
        S.video.writerObj             = VideoWriter(nameVid);
        S.video.frameRate             = 1/S.config.Ts;
        open(S.video.writerObj);
    end
end

function S = gen_x0(S,x_y_theta)
    % DO NOT FORGET CHECKING FEASIBILITY OF INITIAL CONDITION!!!
    if nargin == 1
        Dx                  = S.path.coordinates(1,2)-S.path.coordinates(1,1);
        Dy                  = S.path.coordinates(2,2)-S.path.coordinates(2,1);
        theta0              = atan2(Dy,Dx);
        thetas              = repmat(theta0,S.config.Nt+1,1);
        xy_0                = S.path.coordinates(:,1);
    else        
%         lenXYTh = length(x_y_theta);
%         switch lenXYTh
%             case 3 % x0, y0, theta0, theta_i = theta0
%                 xAux        = x_y_theta(1);
%                 yAux        = x_y_theta(2);
%                 thetas      = repmat(x_y_theta(3),S.config.Nt+1,1);
%                 xy_0        = [xAux;yAux];
%             case 4 % x0, y0, theta0, theta1, theta_i = theta1
%                 xAux        = x_y_theta(1);
%                 yAux        = x_y_theta(2);
                thetas      = [x_y_theta(3:end)'; repmat(x_y_theta(end),S.config.Nt+1-length(x_y_theta(3:end)),1)];
                xy_0        = x_y_theta(1:2)';%[xAux;yAux];
%             otherwise
%                 
%         end        
    end
%     thetas                  = theta0;%repmat(theta0,S.config.Nt+1,1);
%     if rndTheta
%         for i=1:S.config.Nt
%             thetas = [thetas;thetas(end)+(2*rand-1)*pi/2];
% %             thetas(2:end) = thetas(2:end) + (rand(S.config.Nt,1)*4*pi-2*pi) .* (2*S.config.initUncertainty);
%         end
%     else
%         thetas      = repmat(theta0,S.config.Nt+1,1);
%     end
    betas                   = -diff(thetas);    
    %
    if S.config.Nt == 10
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_9                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
        xy_N                = xy_9 - [S.system.Lh10*cos(thetas(10))+S.system.L10*cos(thetas(11)); S.system.Lh10*sin(thetas(10))+S.system.L10*sin(thetas(11))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_8;xy_9;xy_N];
    elseif S.config.Nt == 9
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_N                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_8;xy_N];
    elseif S.config.Nt == 8
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_N                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_N];
    elseif S.config.Nt == 7
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_N                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_N];
    elseif S.config.Nt == 6
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_N                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_N];
    elseif S.config.Nt == 5
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_N                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_N];
    elseif S.config.Nt == 4
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_N                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_N];
    elseif S.config.Nt == 3
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_N                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_N];
    elseif S.config.Nt == 2
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_N                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_0toN = [xy_0;xy_1;xy_N];
    else
        xy_N                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_0toN = [xy_0;xy_N];
    end
    %
    S.init_condition.x0 = [ betas; thetas; xy_0toN];%; zeros(S.system.nu,1) ];
    %
    S.data.xsim(:,1)    = S.init_condition.x0;%[S.init_condition.x0; [1;1]];
%     S.data.xsim(:,1)   = S.path.ref(:,1);%[S.init_condition.x0; [1;1]];
end

function S = gen_x0bar(S)
    % DO NOT FORGET CHECKING FEASIBILITY OF INITIAL CONDITION!!!
    thetas  = repmat(S.init_condition.x0(S.config.Nt+1),S.config.Nt+1,1);
    thetas(2:end) = thetas(2:end) + (rand(S.config.Nt,1)-0.5) .* (2*S.config.initUncertainty);
    betas   = -diff(thetas);
    xy_0    = S.init_condition.x0(2*S.config.Nt+2:2*S.config.Nt+3) + S.config.noise_lvl(S.config.Nt+2:S.config.Nt+3).*randn(2,1);
    %
    if S.config.Nt == 10
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_9                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
        xy_N                = xy_9 - [S.system.Lh10*cos(thetas(10))+S.system.L10*cos(thetas(11)); S.system.Lh10*sin(thetas(10))+S.system.L10*sin(thetas(11))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_8;xy_9;xy_N];
    elseif S.config.Nt == 9
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_N                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_8;xy_N];
    elseif S.config.Nt == 8
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_N                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_N];
    elseif S.config.Nt == 7
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_N                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_N];
    elseif S.config.Nt == 6
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_N                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_N];
    elseif S.config.Nt == 5
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_N                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_N];
    elseif S.config.Nt == 4
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_N                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_N];
    elseif S.config.Nt == 3
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_N                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_N];
    elseif S.config.Nt == 2
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_N                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_0toN = [xy_0;xy_1;xy_N];
    else
        xy_N                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_0toN = [xy_0;xy_N];
    end   
    
    if S.config.SIM
        S.init_condition.x0bar = [ betas; thetas; xy_0toN];%; zeros(S.system.nu,1) ];
% S.init_condition.x0bar = S.init_condition.x0;
    else    
        S                       = read_sensors(S);
        S                       = measurements_vector(S);
        
        xy_0                    = [S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0];
        xy_1                    = xy_0 - [S.system.Lh1*cos(S.ROS.sensors.vectornav_theta0)+S.system.L1*cos(S.ROS.sensors.vectornav_theta0); S.system.Lh1*sin(S.ROS.sensors.vectornav_theta0)+S.system.L1*sin(S.ROS.sensors.vectornav_theta0)];
        xy_N                    = xy_1 - [S.system.Lh2*cos(S.ROS.sensors.vectornav_theta0)+S.system.L2*cos(S.ROS.sensors.vectornav_theta0); S.system.Lh2*sin(S.ROS.sensors.vectornav_theta0)+S.system.L2*sin(S.ROS.sensors.vectornav_theta0)];

        S.init_condition.x0bar = [ S.ROS.sensors.encoders.beta1;    % The vehicle must start aligned, i.e., all segment with the same attitude
                                   S.ROS.sensors.encoders.beta2;
                                   repmat(S.ROS.sensors.vectornav_theta0,S.config.Nt+1,1);...
                                   xy_0;...
                                   xy_1;...
                                   xy_N];
        
    end    
% ####    
S.data.xest = S.init_condition.x0bar;
% ####
end
