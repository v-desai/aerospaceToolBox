%% plotPlanetTEST
close all; clear all; clc

position = [6990077.798814194;1617465.311978378;22679.810569245355];
velocity = [-1675.13972506056;7273.72441330686;252.688512916741];
RV = [position;velocity]*1e-3;
stationAtoll = [-6143584;1364250;1033743]*1e-3;
stationDiego = [1907295;6030810;-817119]*1e-3;
stationArecibo = [2390310;-5564341;1994578]*1e-3;

% a = 6793;
% mu = 398600.4415;
% r = norm([x y z]');
% T = 2*pi*sqrt(r^3/mu);
t0 = 0;
tf = 86400;

input.runTwoBody = 1;
options = odeset('reltol',1e-12,'abstol',1e-12);
[T,RVset] = ode45(@rvProp,[t0 tf],RV ,options, input);
RVset = RVset(:,1:6);

planetName = 'Earth';
solarSystem = 0;
grid = 0;
h1 = plotPlanet(planetName,solarSystem,grid); view(24,90);
laserPlot = [];
hold on;

originalPlot = plot3(RVset(:,1),RVset(:,2),RVset(:,3),'w');
hold on

for ii = 1:10:length(RVset)
    
    tail = 1;
    if ii > 100
        tail = 1+ii-100;
    end
    
    tailPlot = plot3(RVset(tail:ii,1),RVset(tail:ii,2),RVset(tail:ii,3),'b','LineWidth',5);
    hold on
    sc = scatter3(RVset(ii,1),RVset(ii,2),RVset(ii,3),'filled','w');
    hold on
    stationPlot = scatter3(stationArecibo(1),stationArecibo(2),stationArecibo(3),15,'filled','g');
    hold on
    
%     u = stationArecibo;
%     v = RVset(ii,1:3)';
%     ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));
    u = stationArecibo;
    v = RVset(ii,1:3)';
    LOS = v-u;
    
    ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));
    
    if ThetaInDegrees < 30 && ThetaInDegrees > 0
        laserPlot = plot3([RVset(ii,1),stationArecibo(1)],[RVset(ii,2),stationArecibo(2)],[RVset(ii,3),stationArecibo(3)],'b');
        ThetaInDegrees
        T(ii)
    end
    
    hold on
    
    drawnow
    pause(0.001)
    delete(sc)
    delete(tailPlot)
    if ~isempty(laserPlot)
        %delete(laserPlot)
    end
%     view([0-20*(ii/1000),10])
    axis vis3d
    hold on
    
end