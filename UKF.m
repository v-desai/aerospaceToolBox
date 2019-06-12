%% ========== DEMONSTRATION OF BASIC UNSCENTED KALMAN FILTER ==========
% Date:   3/8/2019
% Author: Vivek Desai
% Algorithm verified with Dr. Jones - Orbital Debris: HW 2 Solutions
clear all ; clc ; close all ;
% ------------------------------------------------------------------------

%% ========== CONFIGURATION ==========

% ----- Retrieve Data -----

load testData_angles.dat                % columns [time] [RA] [DEC] [x,y,z] [vx,vy,vz]
time = testData_angles(:,1);            % [seconds] [161,1]
measurements = testData_angles(:,2:3);  % [radians] [161,2]
truthPosVel = testData_angles(:,4:9);   % [km , km/s] [161,6]
measurementType = 'RADEC';

% ----- Set Initial Conditions -----

% A priori
xBar0 = [-2011.990;-382.065;6316.376;5.419783;-5.945319;1.37398];
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

%% ========== RUN UNSCENTED KALMAN FILTER ==========

for kk = 1:length(time)
    
    % ====== Time Update ======
    
    if time(kk) == 0        
        xBar = xBar0;
        PBar = PBar0;
        [sigmaPoints] = UKFgenerateSigmaPoints(xBar,PBar,alpha,k);
    else
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
    
    % ----- Calculate post-fit measurement residuals -----
    
    [sigmaPoints] = UKFgenerateSigmaPoints(xHat,P,alpha,k);
    [sigmaPointsZ] = UKFgenerateSigmaPointsZ(sigmaPoints,measurementType);
    [zBar] = UKFmeanMeasurementEstimate(W0_m,Wi_m,sigmaPointsZ);
    [PBarzz] = UKFmeanMeasurementCovariance(W0_c,Wi_c,zBar,sigmaPointsZ);
    
    % ----- Store Post-fit values -----
    
    store_xHat(:,kk) = xHat;
    store_P(:,:,kk) = P;
    store_zResiduals(:,kk) = z - zBar;
    store_Pbarzz(:,:,kk) = PBarzz;
    
    progress = kk/length(time);
    
    fprintf('\b\b')
    fprintf('%.0f',progress*100)
end

%% ==================== Post-filter Analysis ====================

% ----- Figures -----

timePlot = time./3600;
idxConverge = 1;

plotRADEC(idxConverge,store_prezResiduals,timePlot,'Pre-fit Residuals (\alpha,\delta)')
plotRADEC(idxConverge,store_zResiduals,timePlot,'Post-fit Residuals (\alpha,\delta)')

residuals = store_xHat(1:6,:) - truthPosVel';
sigmaPosVel = zeros(length(time),n);
for ii = 1:length(time)
    sigmaPosVel(ii,1:6) = diag(store_P(1:6,1:6,ii)).^0.5;
end
plotPosVel(idxConverge,residuals',sigmaPosVel,timePlot,'Post-fit Residuals (Cartesian)')

save('saveUKF.mat','store_xBar','store_xHat','store_PBar','store_P','truthPosVel')