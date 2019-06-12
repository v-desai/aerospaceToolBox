%% Solve Legrange Planetary Equations
function dOE = legrangeEquations(OE,mu,dRdOE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Traditional set of Keplerian Orbit Elements (OE)
% - OE(1)       Semi-major axis [m]
% - OE(2)       Eccentricity [no unit]
% - OE(3)       Inclination [degrees]
% - OE(4)       Argument of periapsis [degrees]
% - OE(5)       Longitude of ascending node [degrees]
% - OE(6)       Mean anomaly at t0
% - mu          Standard gravitational parameter of central body [m^3/s^2]
% - dRdOE [6x1] Partials of the disturbing potential w.r.t. OE 
%
% OUTPUTS
% - dOE   [6x1] Time rate of change of OE
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = OE(1);
e = OE(2);
i = OE(3)*pi/180;
omega = OE(4)*pi/180;
Omega = OE(5)*pi/180;
M0 = OE(6)*pi/180;

dRda = dRdOE(1);
dRde = dRdOE(2);
dRdi = dRdOE(3);
dRdomega = dRdOE(4);
dRdOmega = dRdOE(5);
dRdM0 = dRdOE(6);

n = sqrt(mu/a^3);
dadt = 2/(n*a) * dRdM0;
dedt = (1-e^2)/(n*a^2*e)*dRdM0 - sqrt(1-e^2)/(n*a^2*e)*dRdomega;
didt = 1/(n*a^2*sqrt(1-e^2)*sin(i)) * (cos(i)*dRdomega - dRdOmega);
domegadt = sqrt(1-e^2)/(n*a^2*e)*dRde - cot(i)/(n*a^2*sqrt(1-e^2))*dRdi;
dOmegadt = 1/(n*a^2*sqrt(1-e^2)*sin(i))*dRdi;
dM0dt = -(1-e^2)/(n*a^2*e)*dRde - 2/(n*a)*dRda;
dMdt = n + dM0dt;

dOE = [dadt dedt didt*180/pi domegadt*180/pi dOmegadt*180/pi dMdt*180/pi]';