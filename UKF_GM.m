%% ========== DEMONSTRATION OF GM UNSCENTED KALMAN FILTER ==========
% Date:   3/8/2019
% Author: Vivek Desai
% Algorithm verified with Dr. Jones - Orbital Debris: HW 3 Solutions
clear all ; clc ; close all ;
% ------------------------------------------------------------------------

%% ========== CONFIGURATION ==========

% ----- Retrieve Data -----

load testData_angles_GM.dat                % columns [time] [RA] [DEC] [x,y,z] [vx,vy,vz]
time = testData_angles_GM(:,1);            % [seconds] [161,1]
measurements = testData_angles_GM(:,2:3);  % [radians] [161,2]
truthPosVel = testData_angles_GM(:,4:9);   % [km , km/s] [161,6]
measurementType = 'RADEC';
load testData_GMcomponents.mat             % variables [weights] [means] [covars]
componentMeans = means;
componentWeights = weights;
componentCovariances = covars;

% ----- Set Initial Conditions -----

% A priori
xBar0 = [-2011.990 ; -382.065 ; 6316.376 ; 5.419783 ; -5.945319 ; 1.37398];
PBar0 = [eye(3),zeros(3) ; zeros(3),10^-6*eye(3)];
% Process Noise
sigmaQ = 10^-8;                         % [km/s^2]
Q = diag([sigmaQ^2,sigmaQ^2,sigmaQ^2]);
% Measurement Noise
sigmaR = 10*pi/(180*3600);              % [radians]
R = diag([sigmaR^2,sigmaR^2]);
% Save common array lengths
n = size(xBar0,1);
d = size(Q,1);
m = size(R,1);
N = n+d;
% Include Process Noise in State
xBar0 = [xBar0;zeros(3,1)];
PBar0(n+1:N,n+1:N) = Q;
% Propagation
input   = [1 1 1 1 0 1]; % [TwoBody,J2,J3,Drag,STM,UKF];
options = odeset('reltol',1e-12,'abstol',1e-12);
% Sigmapoint constants
Beta = 2; alpha = 1; k = 3-N;
[W0_m,W0_c,Wi_m,Wi_c] = UKFgenerateWeights(N,Beta,alpha,k);
% Gaussian Mixture Model constants
mergingThreshhold = 3;
pruningThreshhold = 10^-5;
maxComponents     = 800;
J = size(componentMeans,2);

% ----- Pre-allocate variables used in filter -----

sigmaPoints = zeros(N,2*N+1);
store_prezResiduals = zeros(m,length(time));
store_prePbarzz = zeros(m,m,length(time));
store_xBar = zeros(N,length(time));
store_PBar = zeros(N,N,length(time));
store_xHat = zeros(N,length(time));
store_P = zeros(N,N,length(time));
store_zResiduals = zeros(m,length(time));
store_Pbarzz = zeros(m,m,length(time));
store_GMcomponents = zeros(1,length(time));
xHatGM = zeros(N,J);
PGM = zeros(N,N,J);
zbarGM = zeros(m,J);
zGM = zeros(m,J);
PBarzzGM = zeros(m,m,J);

%% ========== RUN UNSCENTED KALMAN FILTER ==========

f = waitbar(0,'1','Name','Matlab Application');

for kk = 1:length(time)
       
    for ii = 1:J
        
        % ====== Time Update ======
        
        if time(kk) == 0
            
            xBar = [componentMeans(:,ii);zeros(3,1)];
            PBar = componentCovariances(:,:,ii);
            PBar(7:9,7:9) = Q;
            [W0_m,W0_c,Wi_m,Wi_c] = UKFgenerateWeights(N,Beta,alpha,k);
            [sigmaPoints] = UKFgenerateSigmaPoints(xBar,PBar,alpha,k);
            
        else
            
            [sigmaPoints] = UKFgenerateSigmaPoints(xHatGM(:,ii),PGM(:,:,ii),alpha,k);                        
            for jj = 1:2*N+1
                [T,RVset] = ode45(@rvProp,[time(kk-1) time(kk)],sigmaPoints(:,jj),options,input);
                sigmaPoints(:,jj) = RVset(end,:)';
            end
            [xBar] = UKFmeanStateEstimate(W0_m,Wi_m,sigmaPoints);
            [PBar] = UKFmeanStateCovariance(W0_c,Wi_c,xBar,sigmaPoints);
            
        end
        
        % ====== Measurement Update ======
        
        [sigmaPointsZ] = UKFgenerateSigmaPointsZ(sigmaPoints,measurementType);
        [zBar] = UKFmeanMeasurementEstimate(W0_m,Wi_m,sigmaPointsZ);
        [PBarzz] = UKFmeanMeasurementCovariance(W0_c,Wi_c,zBar,sigmaPointsZ);
        PBarzz = PBarzz + R;
        [PBarxz] = UKFmeanCrossCovarianceMatrix(W0_c,Wi_c,sigmaPoints,xBar,sigmaPointsZ,zBar);
        z = measurements(kk,:)';
        
        % ----- Store pre-fit values -----
        
        store_xBar(:,kk) = xBar;
        store_PBar(:,:,kk) = PBar;
        store_prezResiduals(:,kk) = z - zBar;
        store_prePbarzz(:,:,kk) = PBarzz;
        
        % ----- Update state and covariance -----
        
        [xHat,P] = UKFkalmanUpdate(PBarxz,PBarzz,xBar,PBar,z,zBar,Q,n,d,N);
        
        % ----- Store component values -----
        
        xHatGM(:,ii) = xHat;
        PGM(:,:,ii) = P;
        zbarGM(:,ii) = zBar;
        zGM(:,ii) = z;
        PBarzzGM(:,:,ii) = PBarzz;
        
    end
    
    % ----- Store Pre-fit values -----
    
    [xBar,PBar] = GMupdateMeanCovariance(xHatGM,PGM,componentWeights);
    store_xBar(:,kk) = xBar;
    store_PBar(:,:,kk) = PBar;
    
    % ----- Update weights -----
    
    [wHat] = GMupdateWeights(componentWeights,zGM(:,1:J),zbarGM(:,1:J),PBarzzGM(:,:,1:J));
    
    % ----- Merge/Prune/Truncate Algorithm -----
    
    [componentWeights,xHatGM,PGM] = GMmergeComponents(pruningThreshhold,mergingThreshhold,maxComponents,wHat,xHatGM,PGM);
    J = size(xHatGM,2);

    % ----- Calculte PDF Mean and Covariance -----
    
    [xHat,P] = GMupdateMeanCovariance(xHatGM,PGM,componentWeights);
    
    % ----- Calculate post-fit measurement residuals -----
    
    [sigmaPoints] = UKFgenerateSigmaPoints(xHat,P,alpha,k);
    [sigmaPointsZ] = UKFgenerateSigmaPointsZ(sigmaPoints,'RADEC');
    [zbar] = UKFmeanMeasurementEstimate(W0_m,Wi_m,sigmaPointsZ);
    [Pbarzz] = UKFmeanMeasurementCovariance(W0_c,Wi_c,zbar,sigmaPointsZ);
    
    % ----- Store Post-fit values -----
    
    store_xHat(:,kk) = xHat;
    store_P(:,:,kk) = P;
    store_zResiduals(:,kk) = z - zbar;
    store_Pbarzz(:,:,kk) = Pbarzz;
    store_GMcomponents(:,kk) = J;

    % -------------------------------------------------------------------
    progress = kk/length(time);
    percentile = [sprintf('%.0f%',progress*100),'%'];
    waitbar(progress,f,percentile);
    
end

delete(f)

%% ==================== Post-filter Analysis ====================

% ----- Figures -----

timePlot = time./3600;
idxConverge = 1;

plotRADEC(idxConverge,store_zResiduals,timePlot,'Post-fit Residuals (\alpha,\delta)')

preresiduals = store_xHat(1:6,:) - truthPosVel';
sigmaPosVel = zeros(length(time),n);
for ii = 1:length(time)
    sigmaPosVel(ii,1:6) = diag(store_P(1:6,1:6,ii)).^0.5;
end
plotPosVel(idxConverge,preresiduals',sigmaPosVel,timePlot,'Post-fit Residuals (Cartesian)')

figure()
plot(timePlot,store_GMcomponents,'k')
set(gca, 'YScale', 'log')
xlabel('Time since Epoch [hrs]')
ylabel('Number of Components')
grid on;

save('saveUKF_GM.mat','store_xBar','store_xHat','store_PBar','store_P','truthPosVel')