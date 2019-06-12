%% CW Equations Analytical
function [CW_t0_tk,phi_rr,phi_rv,phi_vr,phi_vv] = cwEquations(n,t)

phi_rr = [4-3*cos(n*t) 0 0 ; 6*(sin(n*t)-n*t) 1 0 ; 0 0 cos(n*t)];
phi_rv = [1/n*sin(n*t) 2/n*(1-cos(n*t)) 0 ; 2/n*(cos(n*t)-1) 1/n*(4*sin(n*t)-3*n*t) 0 ; 0 0 1/n*sin(n*t)];
phi_vr = [3*n*sin(n*t) 0 0 ; 6*n*(cos(n*t)-1) 0 0 ; 0 0 -n*sin(n*t)];
phi_vv = [cos(n*t) 2*sin(n*t) 0 ; -2*sin(n*t) 4*cos(n*t)-3 0 ; 0 0 cos(n*t)];

CW_t0_tk = [ phi_rr , phi_rv ; phi_vr , phi_vv ];