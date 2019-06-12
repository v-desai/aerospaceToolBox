function Beta = IJK2RIC(rvec,vvec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% - rvec [3x1] Position vector in Cartesian [LU]
% - vvec [3x1] Velocity vector in Cartesian [LU/TU]
%
% OUTPUT
% - Beta [3x3] Transformation matrix from Cartesian to Radial, In-track,
%              Crosstrack
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = rvec/norm(rvec);
C = cross(rvec,vvec)/norm(cross(rvec,vvec));
I = cross(C,R);
Beta = [R I C]';