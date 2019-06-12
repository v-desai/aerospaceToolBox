function dRV = propThreeBody(t,X,MU,NL)
%% Circular Restricted Three Body Problem
if NL(1,1) == 1 % Nonlinear
    x = X(1);
    y = X(2);
    z = X(3);
    xdot = X(4);
    ydot = X(5);
    zdot = X(6);
    dRV = EQa_CRTPB_NL(MU,x,xdot,y,ydot,z,zdot);
else % Linearized (given deviations)
    dX = X(1:6,1);
    r0 = NL(2:4,1);
    A = EQA_CRTBP(MU,r0(1),r0(2),r0(3));
    dRV = A*dX;
end


end



