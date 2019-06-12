%% Cartesian State Vectors --> Keplerian Orbit Elements
function [OE] = RV2OE(RV,dt,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Cartesian position and velocity (R,V)
% - RV(6x1) [position vector [m] velocity vector [m/s]]
% - mu    Standard gravitational parameter of central body [m^3/s^2]
%
% OUTPUTS
% Traditional set of Keplerian Orbit Elements (OE)
% - OE(1) Semi-major axis [m]
% - OE(2) Eccentricity [no unit]
% - OE(3) Inclination [degrees]
% - OE(4) Argument of periapsis [degrees]
% - OE(5) Longitude of ascending node [degrees]
% - OE(6) Mean anomaly at t0
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE NAME ASSIGNMENT
rvec = RV(1:3,1);
vvec = RV(4:6,1);
%% Algorithm
% Calculate orbital momentum
hvec = cross(rvec,vvec);
% Calculate eccentricity vector
evec = cross(vvec,hvec)/mu - rvec/norm(rvec);
e = norm(evec);
% Calculate vector pointing towards the ascending node
nvec = cross([0 0 1],hvec);
% Calculate true anomaly
if dot(rvec,vvec) >= 0
    nu = acos(dot(evec,rvec)/(norm(evec)*norm(rvec)));
else
    nu = 2*pi- acos(dot(evec,rvec)/(norm(evec)*norm(rvec)));
end
% Calculate orbit inclination
i = acos(hvec(3)/norm(hvec));
% Calculate eccentric anomaly
E = 2*atan(tan(nu/2)/sqrt((1+e)/(1-e)));
% Calculate longitude of ascending node
if nvec(2) >= 0
    Omega = acos(nvec(1)/norm(nvec));
elseif nvec(2) < 0
    Omega = 2*pi - acos(nvec(1)/norm(nvec));
end
% Calculate argument of periapsis
if evec(3) >= 0
    omega = acos(dot(nvec,evec)/(norm(nvec)*norm(evec)));
elseif evec(3) < 0
    omega = 2*pi - acos(dot(nvec,evec)/(norm(nvec)*norm(evec)));
end
% Compute semi-major axis
a = 1/(2/norm(rvec)-norm(vvec)^2/mu);
% Compute mean anomaly
M0 = E - e*sin(E) - sqrt(mu/a^3)*dt;
while M0 > 2*pi
    M0 = M0 - 2*pi;
end
while M0 < 0
    M0 = M0 + 2*pi;
end
M = E - e*sin(E);

% VARIABLE REASSIGNMENT
OE = [a e i*180/pi omega*180/pi Omega*180/pi M*180/pi]';

% test
% 
% [OE] = CART2OE(RV,dt,mu)
