%% Solve Guass Equations
function dOE = guassEquations(OE_M,STW,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Traditional set of Keplerian Orbit Elements (OE)
% - OE_M(1)       Semi-major axis [m]
% - OE_M(2)       Eccentricity [no unit]
% - OE_M(3)       Inclination [degrees]
% - OE_M(4)       Argument of periapsis [degrees]
% - OE_M(5)       Longitude of ascending node [degrees]
% - OE_M(6)       Mean anomaly at t0
% - mu            Standard gravitational parameter of central body [m^3/s^2]
% - STW           Disturbing acceleration in STW direction
%
% OUTPUTS
% - dOE     [6x1] Time rate of change of OE
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = STW(1);
T = STW(2);
W = STW(3);

a = OE_M(1);
e = OE_M(2);
i = OE_M(3)*pi/180;
omega = OE_M(4)*pi/180;
Omega = OE_M(5)*pi/180;
M = OE_M(6)*pi/180;

E = KeplerSolver(M,e);
nu = trueAnomaly(E,e);

p = a*(1-e^2);
h = sqrt(mu*p);
n = sqrt(mu/a^3);
r = a*(1-e^2)/(1+e*cos(nu));

adot = 2/(n*sqrt(1-e^2))*((e*sin(nu))*S+(1+e*cos(nu))*T);
edot = h/mu*(sin(nu)*S+((e+cos(nu))/(1+e*cos(nu))+cos(nu))*T);
idot = h/mu*(cos(nu+omega))/(1+e*cos(nu))*W;
Omegadot = h/(mu*sin(i))*(sin(nu+omega)/(1+e*cos(nu))*W);
omegadot = h/(mu*e)*((2+e*cos(nu))*sin(nu)/(1+e*cos(nu))*T-cos(nu)*S)-cos(i)*Omegadot;
sigmadot = r/(n*a^2*e)*((cos(nu)+e*cos(nu)^2-2*e)*S-sin(nu)*(2+e*cos(nu))*T);
Mdot = sqrt(mu/a^3) + sigmadot;

dOE = [adot edot idot*180/pi omegadot*180/pi Omegadot*180/pi Mdot*180/pi]';