%% Keplerian Orbit Elements --> Cartesian State Vectors
function [RV] = OE2RV(OE,dt,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Traditional set of Keplerian Orbit Elements (OE)
% - OE(1) Semi-major axis [m]
% - OE(2) Eccentricity [no unit]
% - OE(3) Inclination [degrees]
% - OE(4) Argument of periapsis [degrees]
% - OE(5) Longitude of ascending node [degrees]
% - OE(6) Mean anomaly at t0
% - dt     Time of flight of object (t-t0) [sec]
% - mu    Standard gravitational parameter of central body [m^3/s^2]
%
% OUTPUTS
% Cartesian position and velocity (R,V)
% - R[3x1] Position vector [m]
% - V[3x1] Velocity vector [m/s]
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE NAME ASSIGNMENT
a = OE(1);
e = OE(2);
i = OE(3)*pi/180;
omega = OE(4)*pi/180;
Omega = OE(5)*pi/180;
M0 = OE(6)*pi/180;
%% Algorithm
M = M0 + dt*sqrt(mu/a^3);
E = KeplerSolver(M,e);
% true anomaly
nu = trueAnomaly(E,e);

rabs = a*(1-e^2)/(1+e*cos(nu));
h = (mu*a*(1-e^2))^(1/2);
p=a*(1-e^2);

r(1) = rabs*(cos(Omega)*cos(omega+nu) - sin(Omega)*sin(omega+nu)*cos(i));
r(2) = rabs*(sin(Omega)*cos(omega+nu) + cos(Omega)*sin(omega+nu)*cos(i));
r(3) = rabs*(sin(i)*sin(omega+nu));
v(1) = (r(1)*h*e)/(rabs*p)*sin(nu) - h/rabs*(cos(Omega)*sin(omega+nu) + sin(Omega)*cos(omega+nu)*cos(i));
v(2) = (r(2)*h*e)/(rabs*p)*sin(nu) - h/rabs*(sin(Omega)*sin(omega+nu) - cos(Omega)*cos(omega+nu)*cos(i));
v(3) = (r(3)*h*e)/(rabs*p)*sin(nu) + h/rabs*sin(i)*cos(omega+nu);

RV = [r v]';

% CHECK
% [RV] = OE2RV([1.5,0.6,70,-30,160,0]',0,1)