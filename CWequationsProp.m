%% ODE45 - Propagates the C-W Equations
function dX = CWequationsProp(t,X,mu,R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% - t      Current time [TU]
% - mu     Gravitational parameter of central body [LU^3/TU^2]
% - R      Chief distance from center of central body [LU]
% - X(1)   Deputy position in radial direction [LU]
% - X(2)   Deputy position in in-track direction [LU]
% - X(3)   Deputy position in crosstrack direction [LU]
% - X(4)   Deputy velocity in radial direction [LU/TU]
% - X(5)   Deputy velocity in in-track direction [LU/TU]
% - X(6)   Deputy velocity in crosstrack direction [LU/TU]
%
% OUTPUTS
% - dX(1)  Deputy velocity in radial direction [LU/TU]
% - dX(2)  Deputy velocity in in-track direction [LU/TU]
% - dX(3)  Deputy velocity in crosstrack direction [LU/TU] 
% - dX(4)  Deputy acceleration in radial direction [LU/TU]
% - dX(5)  Deputy acceleration in in-track direction [LU/TU]
% - dX(6)  Deputy acceleration in crosstrack direction [LU/TU] 
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = sqrt(mu/R^3);

x = X(1); xdot = X(4);
y = X(2); ydot = X(5);
z = X(3); zdot = X(6);

dX(1,1) = xdot;
dX(2,1) = ydot;
dX(3,1) = zdot;

dX(4,1) = 3*n^2*x + 2*n*ydot;
dX(5,1) = -2*n*xdot;
dX(6,1) = -n^2*z;