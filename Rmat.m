%% Simple euler rotations
function R = Rmat(axis,angle,oneifrad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% - axis       Rotation axis; x-axis (1) y-axis (2) z-axis (3)
% - angle      Angle of rotatation [radians or degrees]
% - oneifrad   Designates if input angle is in radians (1) or degrees 
%              (other)
%
% OUTPUTS
% - R [3x3]   Rotation matrix of given axis and angle
%
% Author - Vivek Desai
%
% Last Update - 11/20/2018
%
% References 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if oneifrad == 1
    % angle in radians
else
    % angle in degrees
    angle = deg2rad(angle);
end

if axis == 1
    R = [1 0 0 ; 0 cos(angle) sin(angle) ; 0 -sin(angle) cos(angle)];
end
if axis == 2
    R = [cos(angle) 0 -sin(angle) ; 0 1 0 ; sin(angle) 0 cos(angle)];
end
if axis == 3
    R = [cos(angle) sin(angle) 0 ; -sin(angle) cos(angle) 0 ; 0 0 1];
end

