function dRV = twoBody6stateProp(t,RV,mu,input)
dRV = zeros(length(RV),1);
x = RV(1);
y = RV(2);
z = RV(3);
dRV(1:3,1) = RV(4:6,1);

if input(1,1) == 1
    ux = input(2);
    uy = input(3);
    uz = input(4);
    T = input(2);
    dRV(4:6,1) = EQa_TwoBody(T,mu,x,y,z);
    % dRV(4:6,1) = J2EOM(J2,R,mu,RV(1),RV(2),RV(3));
    
    STM = reshape(RV(7:42),6,6);
    STMdot = F_TwoBody(T,mu,x,y,z)*STM;
    dRV(7:42) = reshape(STMdot,36,1);
else
    T = 0;
    dRV(4:6,1) = EQa_TwoBody(T,mu,x,y,z);
end



