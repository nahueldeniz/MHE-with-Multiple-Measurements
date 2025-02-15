% #########################################################################
% 
% MHE with Multiple Measurements
% Author: Nestor. N. Deniz - 2025
%
% #########################################################################notification icon


function S = N_trailersGPSOutliers()
    clear all; clc;
    import casadi.*
    % ---------------------------------------------------------------------
    S       = init();        
    % Init ROS's things ---------------------------------------------------
    if S.config.dataOnline == true; figure('units','normalized','outerposition',[0 0 1 1]); hold on; end
    if ~S.config.SIM == true; S = ROS(S); end
    for num_sims = 1:S.config.NUM_SIMS
        for num_ampli=1:length(S.config.AMPLITUDES)
            for num_density=1:length(S.config.DENSITIES)
                S.data.num_sims = num_sims;
                S = call_init_functions(S,S.config.DENSITIES(num_density),S.config.AMPLITUDES(num_ampli));
                %
                S.config.iters = 1;
                while (S.config.time(end) < S.config.tf)% && ~(S.path.reach_end_mhempc)
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
                        plot_mono2(S, S.data.xsim(1:S.system.nq,end));       
                        if S.config.video
                            frame = getframe();
                            writeVideo(S.video.writerObj, frame)
                        end
                        S.exec_time.t_plt = [S.exec_time.t_plt, toc];
                    end
                    S.exec_time.t_mis = [S.exec_time.t_mis, toc];
                end
                if ~S.config.SIM
                    S = write_ctrl(S, 0, 0); % Stop the vehiclee            
                end
                if S.config.video
                    close(S.video.writerObj);
                end
                % -----------------------------------------------------------------
                S.data.mhempc.performance.xsim{(num_ampli-1)*length(S.config.DENSITIES) + num_density}       = S.data.xsim;            
                S.data.mhempc.performance.ysim{(num_ampli-1)*length(S.config.DENSITIES) + num_density}       = S.data.ysim;
                S.data.mhempc.performance.xest{(num_ampli-1)*length(S.config.DENSITIES) + num_density}       = S.data.xest;
                S.data.mhempc.performance.ctrl{(num_ampli-1)*length(S.config.DENSITIES) + num_density}       = S.mpc.Controls;
                S.data.mhempc.performance.umeas{(num_ampli-1)*length(S.config.DENSITIES) + num_density}      = S.sensors.velocities;
                S.data.mhempc.performance.mpcRefs{(num_ampli-1)*length(S.config.DENSITIES) + num_density}    = S.data.references;
                S.data.mhempc.performance.exec_times{(num_ampli-1)*length(S.config.DENSITIES) + num_density} = S.exec_time;
                if S.config.save_workspace
                    dirBase = '/home/nahuel/Documents/Nahuel/Algorithms/Obstacle-Avidance N-trailers/Simulations/sims Ts50ms/';
                    nameFile = strcat([dirBase,'obsAvoidanceGNT-Nt',num2str(S.config.Nt),'-Nc',num2str(S.config.Nc),'-Np',num2str(S.config.Np),'-numMovObs',num2str(S.config.maxMovingObs),'-numStatObs',num2str(S.config.maxStaticObs),'-mthd-',S.config.CtrlMethod]);
                    save(nameFile);
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
    % Solver ______________________________________________________________
    S.config.solver         = 'casadi';             % options: 'casadi', 'acado'    
    S.config.CtrlMethod     = 'soft_constraints';   % 'soft_constraints', 'fnmppc'
    S.config.integrator     = 'RK4';                % 'euler'
    S.config.EstMethod      = 'mult_meas';           % 'loss', 'bayesian', 'standard', 'mult_meas';
    % Simulation or field experiment (Husky)_______________________________
    S.config.SIM            = true;
    S.config.stab.compFlg   = false;                % "Smplified design of practically stable MPC schemes", R. Comelli et. al.
    S.config.stab.xMin      = -0.25;
    S.config.stab.xMax      = 0.25;
    S.config.stab.yMin      = -0.25;
    S.config.stab.yMax      = 0.25; 
    S.config.Omega_O        = Polyhedron('lb',[S.config.stab.xMin;S.config.stab.yMin],'ub',[S.config.stab.xMax;S.config.stab.yMax]);
    S.config.Omega_O_updt   = S.config.Omega_O;
    %
    S.config.stab.setNc     = 15;
    %
    S.config.Nt             = 1;
    V                       = [-5 -5; -5 5; 5 5; 5 -5]./(1+log10(S.config.Nt));
    S.config.X_N_Omega_I    = Polyhedron(V);
    S.config.X_N_Omega_I_updt = S.config.X_N_Omega_I;
    S.config.verbose        = false;
    %
    S.config.calcMtxConvCoord = true; if ~S.config.calcMtxConvCoord; warning('WARNING!!! The matrix for correcting x-y coordinates is not being computed...'); end;
    S.config.t0             = 0;
    S.config.Ts             = 0.05;
    S.config.tf             = 3500;    
    S.config.same_seed      = false;
    S.config.EXPORT         = false; if ~S.config.EXPORT; warning('WARNING!!! MHE and MPC were not compiled...'); end
    S.config.iters          = 0;
    S.config.time           = 0;
    S.config.NUM_SIMS       = 1;
    S.config.updtAC         = false;
    S.config.outputs        = [1:S.config.Nt+1,2*S.config.Nt+2:2*S.config.Nt+3];    % 2*S.config.Nt+6:2*S.config.Nt+7];
    S.config.Nc             = 15;                                                   % Lenght of the control horizon
    S.config.Np             = 15;                                                   % prediction horizon
    S.config.Ne             = 15;                                                   % Lenght of the estimation window
    S.config.segmentTosteer = 0;
    % Obstacle configurations _____________________________________________    
    S.config.maxMovingObs   = 0;
    S.config.maxStaticObs   = 0;
    S.config.totNumObs      = S.config.maxStaticObs+S.config.maxMovingObs;
    S.config.obsRand        = false;
    S.config.numStObsSgmnt  = 3;                                                    % number of points to be evaluated in the parametric equation of the line
    S.config.numMvObsSgmnt  = 2;
    S.config.maxDistStatObs = 30;                                                   % Static osbtacles beyond this distance are discarded
    S.config.numMeasEll     = 4;    
    S.config.nroFrDetMObs   = S.config.numMeasEll;
    S.config.thrMovingObs   = 0.002;                                                %0.008;
    S.config.distToUpdEll   = 0.2;
    S.config.accelOneMovObs = true;
    % Estimator and Control algorithms to be executed _____________________
    S.config.obs_strategy   = 'gauss';
    S.config.mhe            = true;
    % Disturbances ________________________________________________________
    S.config.noise          = 'gaussian';                                           % 'gaussian' or 'uniform'
    S.config.noise_lvl      = [(1*pi/180).*ones(S.config.Nt, 1); 0.2*pi/180; [0.025; 0.025]];
    S.config.initUncertainty = 100*pi/180;
    S.config.initXY0Err     = 3;
    S.config.xy0FarFromPath = 1;
    S.config.slip           = [1; 1];                                               % 1= no slippage
    S.config.procDist.type  = 'normal';
    S.config.procDist.amp   = 0.*[(0*pi/180).*ones(S.config.Nt,1); zeros(S.config.Nt+1,1); 0.01.*ones(2,1); 0.01.*ones(2,1)];
    S.config.gpsOutlier     = true;
    S.config.gpsOutlBurstCnt = 0;                                                   % 0: no burst, 1 to Ne number of consecutive outliers
    S.config.AMPLITUDES     = 10;                                                    % The number of elemets of this variable determines the number of simulations to be carried out
    S.config.gpsOutlierAmp  = [];
    S.config.gpsOutlierProbab = [];                                                 % [0, 100]
    S.config.IMUbiasDev     = 5/180;
    S.config.model.uncty    = false;
    S.config.model.dev      = 15;                                                   % porcentual value of uncertainty
    S.config.iterCorrecOn   = 200;
    S.config.CtrlNoise_lvl  = [0; 0];
    S.config.DENSITIES      = 49;                                                   % In %, the number of elemets of this variable determines the number of simulations to be carried out
    % Reference velocities ________________________________________________
    S.config.vNr            = 0;
    S.config.vNTarget       = 1.3;
    S.config.timeToReachVNr = 5;                                                    % seconds
    S.config.deltavN        = S.config.vNTarget/(S.config.timeToReachVNr/S.config.Ts + 2*(S.config.Ne+1));
    S.config.disLastTgt     = 0.75;
    % GPS's relative position to centre of vehicle ________________________
    S.config.gps_x          = -0.3;
    S.config.gps_y          = 0;
    S.config.gps_d          = sqrt(S.config.gps_x^2 + S.config.gps_y^2);
    S.config.gps_fcx        = 0;
    S.config.gps_fcy        = 0;
    S.config.SIMGPS         = false;                                                % use to simulate date from GPS when it does not have signal. Just for indoor testing purposes.
    % Obstacles ___________________________________________________________
    S.config.lidarLimits.X  = 5;    %1.0;                                           % distance to points from the lidar
    S.config.lidarLimits.Y  = 10;   %2.3;   
    S.config.lidarLimits.Z  = 1;
    % boundaries for ignoring obstacles
    S.config.x_max          = 13;
    S.config.x_min          = 2;
    S.config.y_max          = 8;
    S.config.y_min          = 1.5;
    S.config.lidarOn        = false;
    S.config.obsDetection   = false;
    % Multiple measurements configuration _________________________________
    S.config.numSensors     = 10;
    S.config.fractionTs     = 10;                                                   % Value in [0, 99]: if measurements are taken in burst, then they are taken in this interval of Ts
    % Plot data online during experiments. USeful for debbuging purposes __
    S.config.dataOnline     = true;
    S.config.save_workspace = false;    
    S.config.video          = false;
    %    
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
    %
    S.estimator.mtxP = [S.estimator.mtxP, [trace(inv(full(S.mheLoss.nlmheLoss.mtxP))); norm(inv(full(S.mheLoss.nlmheLoss.mtxP)))]];
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
    if strcmp(S.config.CtrlMethod,'soft_constraints')
        if ~isempty(S.mpc.mpcCasadi.Qtraj)
            Qtraj = S.mpc.mpcCasadi.Qtraj;
        else
            Qtraj = repmat(S.mpc.mpcCasadi.q0bar, 1, S.config.Np);
        end
    end
    if strcmp(S.config.CtrlMethod,'fnmppc')
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
    tractorBodyPlt.plot('color',clrTractor);
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
            trailerLoadPlt.plot('color',clrLastTrailerLoad);
        else
            trailerLoadPlt.plot('color',clrTrailerLoad);
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

function S = call_init_functions(S,outlDensity,outlAmplitude)
    S = reserve_temp_memory(S);
    S = gen_init_conditions(S);
    S = init_flags_and_counters(S);
    % Update density and amplitude of the outliers
    S.config.gpsOutlierAmp      = outlAmplitude;
    S.config.gpsOutlierProbab   = outlDensity;
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
            noise = randn(S.system.ny*S.config.numSensors,1).*repmat(S.config.noise_lvl,S.config.numSensors,1);
        else
            noise = (2*rand(S.system.ny*S.config.numSensors,1)-1).*repmat(S.config.noise_lvl,S.config.numSensors,1);
        end
        if S.config.gpsOutlier == true
            for i=1:S.config.numSensors
                if randi(100) <= S.config.gpsOutlierProbab
                    if strcmp(S.config.noise,'gaussian')
                        noise((i-1)*S.system.ny+S.config.Nt+2:(i-1)*S.system.ny+S.config.Nt+3) = S.config.gpsOutlierAmp*randn(2,1);
                    else
                        noise((i-1)*S.system.ny+S.config.Nt+2:(i-1)*S.system.ny+S.config.Nt+3) = S.config.gpsOutlierAmp*(2.*rand(2,1)-1);
                    end
                end
            end
        end            
        S.data.measNoise    = [S.data.measNoise, noise];      
        ysim = [];
        for i=1:S.config.numSensors
            ysimAux = S.data.xsim(S.config.outputs+S.system.nq*(i-1),end);
            ysim    = [ysim; ysimAux];
        end
        S.data.ysim         = [S.data.ysim, ysim + noise];
        %
        unoise              = randn(S.system.nu,1).*S.config.CtrlNoise_lvl;
        S.data.UmeasNoise   = [S.data.UmeasNoise, unoise];
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
    if S.config.mhe
        if strcmp(S.config.EstMethod,'loss')
            updateMeasurement(S.mheLoss.nlmheLoss, double(S.data.ysim(:,end)));    % in S.data.ysim(:,end-2) I save the velocities too
            updateInput(S.mheLoss.nlmheLoss, double(S.sensors.velocities(:,end)));
        elseif strcmp(S.config.EstMethod,'bayesian')
            updateMeasurement(S.mheAle.nlmheAle, double(S.data.ysim(:,end)));
            updateInput(S.mheAle.nlmheAle, double(S.sensors.velocities(:,end)));   
        elseif strcmp(S.config.EstMethod,'standard')
            updateMeasurement(S.algorithms.mheStandard, double(S.data.ysim(:,end)));
            updateInput(S.algorithms.mheStandard, double(S.sensors.velocities(:,end)));   
        elseif strcmp(S.config.EstMethod,'mult_meas')
            updateMeasurement(S.mheLoss.nlmheMultiMeas, double(S.data.ysim(:,end)));    % in S.data.ysim(:,end-2) I save the velocities too
            updateInput(S.mheLoss.nlmheMultiMeas, double(S.sensors.velocities(:,end)));
        end
    else
        updateMeasurement(S.algorithms.ekfCasadi, double(S.data.ysim(1:end-2,end)));
        updateInput(S.algorithms.ekfCasadi, double(S.sensors.velocities(:,end)));
    end
end

function S = call_ESTIMATOR(S)
    tic;
    if S.config.mhe
        if strcmp(S.config.EstMethod,'loss')
            solve(S.mheLoss.nlmheLoss);
            q_k = S.mheLoss.nlmheLoss.qk;
        elseif strcmp(S.config.EstMethod,'bayesian')
            solve(S.mheAle.nlmheAle);
            q_k = S.mheAle.nlmheAle.x_k;
        elseif strcmp(S.config.EstMethod,'standard')
            solve(S.algorithms.mheStandard);
            q_k = S.algorithms.mheStandard.q_k;
        elseif strcmp(S.config.EstMethod,'mult_meas')
            solve(S.mheLoss.nlmheMultiMeas);
            q_k = S.mheLoss.nlmheMultiMeas.qk;
        end
    else
        solve(S.algorithms.ekfCasadi);
        q_k = S.algorithms.ekfCasadi.x_k;
    end
% q_k = S.data.xsim(:,end);
    S.data.xest = [S.data.xest, q_k];
    S.exec_time.t_mhe  = [S.exec_time.t_mhe, toc];
    %

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
    if (S.path.acumIndx+1)>=kMax
        S.path.last_tgt = true;
    end
%     if isempty(i)
%         i = kMax;
%     end
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

function S = call_PFA(S)
    tic;
    % Handle moving target ================================================
    S = updateTarget(S);
S.controller.ref.xReachable = S.controller.ref.x;
S.controller.ref.yReachable = S.controller.ref.y;
S.controller.ref.thetaReachable = S.controller.ref.theta;    
%     % Handle target's velocity ============================================
%     S = handleTargetVelocity(S);
%     % Handle weighting matrices ===========================================
%     S = handleWeightinMatrices(S);
%     % Handle obstacles ====================================================
%     S = handleObstacles(S);
%     % Handle occlusion ====================================================
%     S = handleOcclusion(S);
%     % Update Omega_O & X_N_Omega_I ========================================
%     S = updateSets(S);
%     % Compute reachable target ============================================
%     S = cmptReachableTgt(S);
%     %
    S.exec_time.t_pfa = [S.exec_time.t_pfa, toc];
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

function S = call_SYSTEM(S)
    if S.config.SIM == true
        simulation_input.x = S.data.xsim(1:S.system.nq,end);
        simulation_input.u = S.mpc.Controls(:,end);
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
        states          = FY1(full(simulation_input.x), full(simulation_input.u)) + simulation_input.slip + simulation_input.dist;   % generated mex file
        for i=1:S.config.numSensors-1
            statesAux   = FYNs(full(simulation_input.x), full(simulation_input.u)) + simulation_input.slip + simulation_input.dist;   % generated mex file
            states      = [states; statesAux];
        end
        S.data.xsim = [S.data.xsim, full(states)];
    else
        S = write_ctrl(S, S.mpc.Husky_w0, S.mpc.Husky_v0);
    end
end

function S = call_CONTROLLER(S)
    if S.config.vNr < S.config.vNTarget
        S.config.vNr = S.config.vNr + S.config.deltavN;
    end
    if strcmp(S.config.CtrlMethod,'soft_constraints')
        S = call_MPC(S);
    elseif strcmp(S.config.CtrlMethod,'fnmppc')
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
    % Output of the system ------------------------------------------------
    S.dynamic.h_rhs  = S.dynamic.q(S.config.outputs);
    S.dynamic.h      = casadi.Function('h', {S.dynamic.q}, {S.dynamic.h_rhs});
    % compute jacobians ---------------------------------------------------
    S.dynamic.jac_fx = casadi.Function('Func_A', {S.dynamic.q, S.dynamic.u}, {S.dynamic.f_rhs.jacobian(S.dynamic.q)}, {'x', 'u'}, {'A'});
    S.dynamic.jac_fu = casadi.Function('Func_B', {S.dynamic.q, S.dynamic.u}, {S.dynamic.f_rhs.jacobian(S.dynamic.u)}, {'x', 'u'}, {'B'});
    S.dynamic.jac_hx = casadi.Function('Func_C', {S.dynamic.q}, {S.dynamic.h_rhs.jacobian(S.dynamic.q)}, {'x'}, {'C'});
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
    % Integrator for simulating the multiple measurements _________________
    Ts2             = S.config.Ts*(S.config.fractionTs/(100*max([1,(S.config.numSensors-1)])));
    Ts1             = S.config.Ts*(100-S.config.fractionTs)/100 + Ts2/2;    
    k1              = S.dynamic.f(S.dynamic.q, S.dynamic.u);
    k2              = S.dynamic.f(S.dynamic.q + Ts1 / 2 * k1, S.dynamic.u);
    k3              = S.dynamic.f(S.dynamic.q + Ts1 / 2 * k2, S.dynamic.u);
    k4              = S.dynamic.f(S.dynamic.q + Ts1 * k3, S.dynamic.u);
    x_rk4           = S.dynamic.q + Ts1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    S.dynamic.FY1   = casadi.Function('FNt', {S.dynamic.q, S.dynamic.u}, {x_rk4});
    S.dynamic.FY1.generate('FY1.c',opts);
    mex FY1.c -largeArrayDims;
    %
    k1              = S.dynamic.f(S.dynamic.q, S.dynamic.u);
    k2              = S.dynamic.f(S.dynamic.q + Ts2 / 2 * k1, S.dynamic.u);
    k3              = S.dynamic.f(S.dynamic.q + Ts2 / 2 * k2, S.dynamic.u);
    k4              = S.dynamic.f(S.dynamic.q + Ts2 * k3, S.dynamic.u);
    x_rk4           = S.dynamic.q + Ts2 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    S.dynamic.FYNs  = casadi.Function('FNt', {S.dynamic.q, S.dynamic.u}, {x_rk4});
    S.dynamic.FYNs.generate('FYNs.c',opts);
    mex FYNs.c -largeArrayDims;
end

function S = init_mheStandard(S)
    S.mpc.box_constraints               = struct;
    S.mpc.box_constraints.QluBounds     = [];
    S.mpc.box_constraints.WluBounds     = [-0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1];
    S.mpc.box_constraints.VluBounds     = [-0.1 0.1; -0.1 0.1; -inf inf; -inf inf];
    S.mpc.box_constraints.ZluBounds     = [repmat(-0.01,S.config.Nt,1), repmat(0.01,S.config.Nt,1)];
    S.mpc.box_constraints.UluBounds     = [];
    %
    S.mpc.Mtxs                          = struct;
    
    S.mpc.Mtxs.Q    = (1e6).*eye(7);%diag([0.5.*ones(1,S.config.Nt),1,0.25.*ones(1,S.config.Nt),ones(1,2),0.25.*ones(1,2*S.config.Nt)]);%,0.2,0.2]);
    S.mpc.Mtxs.R    = diag([0.01 0.01 100 100]);%diag([0.25.*ones(1,S.config.Nt),0.65,1,1]);%,1,1]);
    S.mpc.Mtxs.P    = 1e6.*eye(S.system.nq);
    S.mpc.Mtxs.Z    = 1e0.*eye(S.config.Nt);
    S.mpc.Mtxs.dU   = eye(S.system.nu);
    % nlmheCasadiNt(Ne,x,u,Nt,f_rhs,h_rhs,Mtxs,nq,nu,ny,boxConst,q0bar,Ts,dimensions)
    S.algorithms.mheStandard = nlmheCasadiNt(S.config.Ne,S.dynamic.q,S.dynamic.u,S.config.Nt,S.dynamic.f_rhs,S.dynamic.h_rhs,...
        S.mpc.Mtxs,S.system.nq,S.system.nu,S.system.ny,S.mpc.box_constraints,S.init_condition.x0bar,S.config.Ts,S.system.Lh);
    %
    setSigma(S.algorithms.mheStandard, 1);
    setC(S.algorithms.mheStandard, 1e7);

    setJacobians(S.algorithms.mheStandard,S.dynamic.jac_fx,S.dynamic.jac_fu,S.dynamic.jac_hx);        
end

function S = fill_mheStandard(S)
    for imhe=1:S.config.Ne
        updateMeasurement(S.algorithms.mheStandard, full(S.dynamic.h(S.init_condition.x0bar)));
        updateInput(S.algorithms.mheStandard,zeros(S.system.nu,1));       
    end
    updateMeasurement(S.algorithms.mheStandard, full(S.dynamic.h(S.init_condition.x0bar)));
end

function S = init_mhLoss(S)
    dims        = {};
    dims.nq     = S.system.nq;
    dims.nu     = S.system.nu;
    dims.ny     = S.system.ny;
    boxConst    = [];
    % S.mheLoss.nlmheLoss  = mheOptiLoss(S.config.Ne,S.config.Nt,S.dynamic.FNt,S.dynamic.h,dims,boxConst,S.config.Ts);
    S.mheLoss.nlmheLoss  = mheOptiLoss3(S.config.Ne,S.config.Nt,S.dynamic.FNt,S.dynamic.h,dims,boxConst,S.config.Ts);

    set_x0bar(S.mheLoss.nlmheLoss,S.init_condition.x0bar);
    
    % setSigmaP(S.mheLoss.nlmheLoss,1.50); % 0.03
    % setCP(S.mheLoss.nlmheLoss,1e9);

    setJacobians(S.mheLoss.nlmheLoss,S.dynamic.jac_fx,S.dynamic.jac_fu,S.dynamic.jac_hx);
    
    set_cPrm(S.mheLoss.nlmheLoss,0.15);
    set_alpha(S.mheLoss.nlmheLoss,1);
    
    setMtxW(S.mheLoss.nlmheLoss,(1e-6).*eye(7));        % Process Noise Variance Matrix diag([0.5.*ones(1,S.config.Nt),0.1,0.25.*ones(1,S.config.Nt),0.1,0.1,0.1,0.1])
    setMtxR(S.mheLoss.nlmheLoss,diag([0.01 0.01 100 100]));                              % Measurement Noise Variance Matrix diag(diagmx(0.1.*eye(S.config.Nt+1),10.*eye(2)))
    setVehicleDims(S.mheLoss.nlmheLoss,S.system.Lh);
end

function S = fill_mheLoss(S)
    for imhe=1:S.config.Ne
        updateMeasurement(S.mheLoss.nlmheLoss, full(S.dynamic.h(S.init_condition.x0bar)));
        updateInput(S.mheLoss.nlmheLoss,zeros(S.system.nu,1));       
    end
    updateMeasurement(S.mheLoss.nlmheLoss, full(S.dynamic.h(S.init_condition.x0bar)));
end

function S = init_mheMultiMeas(S)
    dims        = {};
    dims.nq     = S.system.nq;
    dims.nu     = S.system.nu;
    dims.ny     = S.system.ny;
    boxConst    = [];
    % S.mheLoss.nlmheLoss  = mheOptiLoss(S.config.Ne,S.config.Nt,S.dynamic.FNt,S.dynamic.h,dims,boxConst,S.config.Ts);
    S.mheLoss.nlmheMultiMeas  = mheMultiMeas(S.config.Ne,S.config.Nt,S.dynamic.FNt,S.dynamic.h,dims,boxConst,S.config.Ts,S.config.numSensors);

    set_x0bar(S.mheLoss.nlmheMultiMeas,S.init_condition.x0bar);
    
    % setSigmaP(S.mheLoss.nlmheLoss,1.50); % 0.03
    % setCP(S.mheLoss.nlmheLoss,1e9);

    setJacobians(S.mheLoss.nlmheMultiMeas,S.dynamic.jac_fx,S.dynamic.jac_fu,S.dynamic.jac_hx);
    
    set_cPrm(S.mheLoss.nlmheMultiMeas,0.15);
    set_alpha(S.mheLoss.nlmheMultiMeas,1);
    
    setMtxW(S.mheLoss.nlmheMultiMeas,(1e-6).*eye(7));        % Process Noise Variance Matrix diag([0.5.*ones(1,S.config.Nt),0.1,0.25.*ones(1,S.config.Nt),0.1,0.1,0.1,0.1])
    setMtxR(S.mheLoss.nlmheMultiMeas,diag([0.01 0.01 100 100]));                              % Measurement Noise Variance Matrix diag(diagmx(0.1.*eye(S.config.Nt+1),10.*eye(2)))
    setVehicleDims(S.mheLoss.nlmheMultiMeas,S.system.Lh);
end

function S = fill_mheMultiMeas(S)
    for imhe=1:S.config.Ne
        updateMeasurement(S.mheLoss.nlmheMultiMeas, repmat(full(S.dynamic.h(S.init_condition.x0bar)),S.config.numSensors,1));
        updateInput(S.mheLoss.nlmheMultiMeas,zeros(S.system.nu,1));       
    end
    updateMeasurement(S.mheLoss.nlmheMultiMeas, repmat(full(S.dynamic.h(S.init_condition.x0bar)),S.config.numSensors,1));
end

function S = init_mheBayesian(S)
    S.mheAle.l_x     = 1e6;
    S.mheAle.P       = S.mheAle.l_x.*casadi.DM.eye(S.system.nq);           % init. val. of the arrival-cost weight matrix
    S.mheAle.Qx      = 15.*eye(S.system.nq);                                   % states stage-cost
    S.mheAle.Rx      = 10.*eye(S.system.ny);                                   % measurement noise stage-cost
    S.mheAle.Tsmhe   = S.config.Ts;                                     % If 0, mhe does not integrate the systems
    S.mheAle.beta    = 1;
    % *************************************************************
    S.mheAle.nlmheAle   = nlmheAlessandri(S.config.Ne,S.config.Nt,S.mheAle.beta,S.dynamic.FNt,S.dynamic.h,...
        S.mheAle.Qx,S.mheAle.Rx,S.system.nq,S.system.nu,S.system.ny,S.system.nq,S.mheAle.P,[],...
        [],[],[],[],S.config.Ts);
    %
    set_x0bar(S.mheAle.nlmheAle,S.init_condition.x0bar);
    setJacobians(S.mheAle.nlmheAle,S.dynamic.jac_fx,S.dynamic.jac_fu,S.dynamic.jac_hx);

    setVehicleDims(S.mheAle.nlmheAle,S.system.Lh);
end

function S = fill_mheBayesian(S)
    for imhe=1:S.config.Ne
        updateMeasurement(S.mheAle.nlmheAle, full(S.dynamic.h(S.init_condition.x0bar)));
        updateInput(S.mheAle.nlmheAle,zeros(S.system.nu,1));       
    end
    updateMeasurement(S.mheAle.nlmheAle, full(S.dynamic.h(S.init_condition.x0bar)));
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
    S.mpc.box_constraints.QNluBounds    = [repmat([-120*pi/180 120*pi/180],S.config.Nt,1); repmat([-inf inf],S.config.Nt+1,1); repmat([-inf inf],2*(S.config.Nt+1),1)];%; [-1 1; -5 5]];

    S.mpc.box_constraints.UluBounds     = [-2 2;...
                                           -2 2];
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
    % Estimators
    S   = init_mhLoss(S);
    S   = fill_mheLoss(S);
    S   = init_mheBayesian(S);
    S   = fill_mheBayesian(S);
    S   = init_mheStandard(S);
    S   = fill_mheStandard(S);
    S   = init_mheMultiMeas(S);
    S   = fill_mheMultiMeas(S);
    % Controller
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
    S.estimator.mtxP                         = [];
    %
    if S.config.video && ~S.config.stab.compFlg
        nameVid = strcat(['obsAvoidanceGNT-Nt',num2str(S.config.Nt),'-Nc',num2str(S.config.Nc),'-Np',num2str(S.config.Np),'-numMovObs',num2str(S.config.maxMovingObs),'-numStatObs',num2str(S.config.maxStaticObs),'-mthd-',S.config.CtrlMethod,'.avi']);
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
        xy_0                = S.path.coordinates(:,1) + S.config.xy0FarFromPath.*randn(2,1);
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
    S.init_condition.x0 = [ betas; thetas; xy_0toN];
    %
    S.data.xsim(:,1)    = repmat(S.init_condition.x0,S.config.numSensors,1);
end

function S = gen_x0bar(S)
    % DO NOT FORGET CHECKING FEASIBILITY OF INITIAL CONDITION!!!
    thetas  = repmat(S.init_condition.x0(S.config.Nt+1),S.config.Nt+1,1);
    thetas(2:end) = thetas(2:end) + (rand(S.config.Nt,1)-0.5) .* (2*S.config.initUncertainty);
    betas   = -diff(thetas);
    xy_0    = S.init_condition.x0(2*S.config.Nt+2:2*S.config.Nt+3) + S.config.noise_lvl(S.config.Nt+2:S.config.Nt+3).*randn(2,1) + S.config.initXY0Err.*randn(2,1);
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
    % Init these values just for avoiding probelms at solving mpc for first
    % time
    S.controller.ref.xReachable     = xy_0(1);
    S.controller.ref.yReachable     = xy_0(2);
    S.controller.ref.thetaReachable = thetas(1);
    % ---------------------------------------------------------------------

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
