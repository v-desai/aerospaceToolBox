%% True Anomaly Equation
function nu = trueAnomaly(E,e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% - E    Eccentric Anomoly [radians]
% - e    Eccentricity
%
% OUTPUTS
% - nu   True anomaly [radians]
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2));