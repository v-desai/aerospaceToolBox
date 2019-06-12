% UKF
clear all; clc; close all;
% ===== TEST CASE =============

runCase = 4;

switch runCase
    case 1
        filename = 'saveUKF';
    case 2
        filename = 'saveUKF_GM';
    case 3
        filename = 'saveUKF_PDA';
    case 4
        filename = 'azel_od_output_1.json';
    case 5
        filename = 'radec_od_output_1.json';
    case 6
        filename = 'radar_od_output_1.json';
end

[xS,PS,xpost,time] = UKFrtsSmoother(filename);

%% ========== Smoother Analysis ==========


%% ========== Plot Results ===============
figure()
residuals = abs(xpost(1:3,:) - xS(1:3,:));
plot(time,residuals)
legend('x','y','z')
xlabel('Time')
ylabel('Residuals (post - smoothed)')

% time = 0:T-1;
% timePlot = time./3600;

% idxConverge = 1;

% if runCase < 4
%     residuals = xpost(1:6,:) - truthPosVel';
%     sigmaPosVel = zeros(length(time),N);
%     for ii = 1:length(time)
%         sigmaPosVel(ii,1:6) = diag(store_P(1:6,1:6,ii)).^0.5;
%     end
%     [RMSpos,RMSvel] = plotPosVel(idxConverge,residuals',sigmaPosVel,timePlot,'Pre-Smoothed Residuals (Cartesian)');
%
%     residuals = xS(1:6,:) - truthPosVel';
%     sigmaPosVel = zeros(length(time),N);
%     for ii = 1:length(time)
%         sigmaPosVel(ii,1:6) = diag(store_P(1:6,1:6,ii)).^0.5;
%     end
%     [RMSposS,RMSvelS] = plotPosVel(idxConverge,residuals',sigmaPosVel,timePlot,'Post-Smoothed Residuals (Cartesian)');
%
%     diffRMSpos = RMSpos-RMSposS;
% end

