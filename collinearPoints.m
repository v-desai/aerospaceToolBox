function f = collinearPoints(x)
%% Find the L1, L2, L3 points of RTBP
MU = 0.01215;
syms x 
f = -x == -(1-MU)*(x+MU)/(sqrt((x+MU)^2))^3 - MU*(x-1+MU)/(sqrt((x-1+MU)^2))^3;
lagrangePoint = double(solve(f,x));
end
