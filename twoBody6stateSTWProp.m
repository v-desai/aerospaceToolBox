function dRV = twoBody6stateSTWProp(t,RV,mu,STW)
dRV(1:3,1) = RV(4:6,1);
dRV(4:6,1) = -RV(1:3,1)*mu/norm(RV(1:3,1))^3;

r = RV(1:3)';
v = RV(4:6)';
R = r/norm(r);
C = cross(r,v)/norm(cross(r,v));
I = cross(C,R);
R_r0h2ijk = [R;I;C]';

a_r0h = [STW(1);STW(2);STW(3)]; 
a_ijk = R_r0h2ijk*a_r0h;
dRV(4:6,1) = dRV(4:6,1) + a_ijk;