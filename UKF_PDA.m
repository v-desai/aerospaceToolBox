%% ========== DEMONSTRATION OF PDA UNSCENTED KALMAN FILTER ==========
% Date:   2/28/2019
% Author: Vivek Desai
clear all ; clc ; close all ;
% ------------------------------------------------------------------------

%% ========== CONFIGURATION ==========

% ----- Retrieve Data -----

load testData_angles_PDA.mat % variables [actualObsNum] [observations] [P_apriori] [times] [truStates] [X_apriori]
time = times;                % [seconds] [161,1]
measurements = observations; % [radians] [2,3,161]
truthPosVel = truStates';    % [km , km/s] [161,6]
measurementType = 'RADEC';

% ----- Set Initial Conditions -----

% A priori
xBar0 = X_apriori;
PBar0 = P_apriori;
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
input.runTwoBody = 1;
input.runJ2 = 1;
input.runJ3 = 1;
input.runDrag = 1;
input.runDrag(2) = 2; % Cd
input.runDrag(3) = 3.6 / (1000^2); % area km^2
input.runDrag(4) = 1350; % mass kg
input.runDrag(5) = 7.29211585530066e-5; % theta dot rad/s
input.runUKF = 1;
options = odeset('reltol',1e-9,'abstol',1e-9);
% Sigmapoint constants
Beta = 2; alpha = 1; k = 3-N;
[W0_m,W0_c,Wi_m,Wi_c] = UKFgenerateWeights(N,Beta,alpha,k);
% PDA constants
Pg = 0.99;      % Probability of gate
Pd = 0.9;       % Probability of detection
lambdac = 2;    % Poisson

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

f = waitbar(0,'1','Name','Matlab Application');

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
    
    % ----- Probabilistic Data Association -----
    
    [p0,pj,zTilde,zTildej,M] = UKFassignmentPDA(measurements,zBar,PBarzz,lambdac,Pd,Pg,m,kk);
    
    % ----- Store pre-fit values -----
    
    store_xBar(:,kk) = xBar;
    store_PBar(:,:,kk) = PBar;
    store_prezResiduals(:,kk) = zTilde;
    store_prePbarzz(:,:,kk) = PBarzz;
    
    % ----- Update state and covariance -----
    
    [xHat,P] = UKFkalmanUpdatePDA(pj,p0,PBarxz,PBarzz,xBar,PBar,zTilde,zTildej,Q,n,d,N,M);
    
    % ----- Calculate post-fit measurement residuals -----
    
    [sigmaPoints] = UKFgenerateSigmaPoints(xHat,P,alpha,k);
    [sigmaPointsZ] = UKFgenerateSigmaPointsZ(sigmaPoints,measurementType);
    [zBar] = UKFmeanMeasurementEstimate(W0_m,Wi_m,sigmaPointsZ);
    [PBarzz] = UKFmeanMeasurementCovariance(W0_c,Wi_c,zBar,sigmaPointsZ);
    zTildej = zeros(m,M);
    zTilde = zeros(m,1);
    for jj = 1:M
        zTildej(:,jj) = measurements(:,jj,kk) - zBar;
    end
    for jj = 1:M
        zTilde = zTilde + pj(jj)*zTildej(:,jj);
    end
    
    % ----- Store Post-fit values -----
    
    store_xHat(:,kk) = xHat;
    store_P(:,:,kk) = P;
    store_zResiduals(:,kk) = zTilde;
    store_Pbarzz(:,:,kk) = PBarzz;
    
    progress = kk/length(time);
    percentile = [sprintf('%.0f%',progress*100),'%'];
    waitbar(progress,f,percentile);
    
end

delete(f)

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

save('saveUKF_PDA.mat','store_xBar','store_xHat','store_PBar','store_P','truthPosVel')