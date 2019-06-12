%% Solves Kepler's equations
function E = KeplerSolver(M,e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% - M   Mean anomaly [radians]
% - e   Eccentricity
%
% OUTPUTS
% - E   Eccentric Anomaly [radians]
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
eps = 1.0;
E = M;
while (eps>1.e-14)
    iter = iter+1;
    Mn = E-e*sin(E);
    Enew = E + (M-Mn)/(1-e*cos(E));
    
    eps = abs((Enew-E)/max(E,1));
    E = Enew;
    if iter>100
        disp('WARNING: Kepler solver iteration divergence')
        pause
    end
end