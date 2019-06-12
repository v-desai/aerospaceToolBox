% Analyze Stuff For Fun

% x0 = [-2012.151 ; -381.450 ; 6316.615 ; 5.400336 ; -5.916814 ; 1.362965];
% 
% t0 = 0;
% tf = 86400;
% dt = 30;
% 
% runTwoBody = 1;
% runJ2 = 1;
% runJ3 = 1;
% runDrag = 0;
% runSTM = 0;
% runUKF = 0;
% runAEGIS = 0;
% input = [runTwoBody,runJ2,runJ3,runDrag,runSTM,runUKF,runAEGIS];
% options = odeset('reltol',1e-13,'abstol',1e-13);
% 
% [~,RV] = ode45(@rvProp,t0:dt:tf,x0,options,input);
% 
% planetPlot('Earth',0,0)
% hold on
% plot3(RV(:,1).*1000,RV(:,2).*1000,RV(:,3).*1000)

mu = 398600.4415
% v = 0.9679359017969791
% r = 1.7962261229376726
% rVec = [1.02300103 1.07600005 1.01100003]'
% vVec = [ 0.62 0.69999996 -0.24999993]'

% hVec = cross(rVec,vVec)
% eVec = cross(vVec,hVec)/mu - rVec/r
% eVec = 1/mu*((v^2-mu/r)*rVec-dot(rVec,vVec)*vVec)

rVec = [6524.84,6862.875,6448.296]
vVec = [4.901327,5.533756,-1.976341]

ER = 6378.1363;
TU = sqrt(ER^3/mu);
rVec = rVec/ER;
vVec = vVec/ER*TU;

hVec = cross(rVec,vVec)
evec = cross(vVec,hVec)/mu - rVec/norm(rVec)

eVec = 1/mu*((v^2-mu/r)*rVec-dot(rVec,vVec)*vVec)