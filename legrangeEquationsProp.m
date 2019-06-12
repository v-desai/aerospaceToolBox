function dOE = legrangeEquationsProp(t,OE,mu,testCase,J2,R)
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
M = OE(6)*pi/180;

switch testCase
    case 0
        % No perturbing potential
        dRdOE = [0 0 0 0 0 0]';
        dOE = legrangeEquations(OE,mu,dRdOE);
    case 1
        % J2 perturbing potential
        dRdOE = dRdOE_Analytical_J2full(J2,M,R,a,e,i,mu,omega);
        dOE = legrangeEquations(OE,mu,dRdOE);
    case 2
        %J2 average perturbing potential
        dRdOE = dRdOE_Analytical_J2avg(J2,R,a,e,i,mu);
        dOE = legrangeEquations(OE,mu,dRdOE);
end
