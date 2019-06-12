function F_TwoBody = EQF_TwoBody(T,mu,x,y,z)
%EQF_TWOBODY
%    F_TWOBODY = EQF_TWOBODY(T,MU,X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    20-Nov-2018 22:53:32

t2 = abs(x);
t3 = abs(y);
t4 = abs(z);
t5 = t2.^2;
t6 = t3.^2;
t7 = t4.^2;
t8 = t5+t6+t7;
t9 = 1.0./t8.^(3.0./2.0);
t10 = sign(x);
t11 = sign(y);
t12 = 1.0./t8.^(5.0./2.0);
t13 = sign(z);
t14 = 1.0./sqrt(t8);
t15 = T.*t14;
F_TwoBody = reshape([0.0,0.0,0.0,t15-mu.*t9-T.*t2.*t9.*t10.*x+mu.*t2.*t10.*t12.*x.*3.0,-T.*t2.*t9.*t10.*y+mu.*t2.*t10.*t12.*y.*3.0,-T.*t2.*t9.*t10.*z+mu.*t2.*t10.*t12.*z.*3.0,0.0,0.0,0.0,-T.*t3.*t9.*t11.*x+mu.*t3.*t11.*t12.*x.*3.0,t15-mu.*t9-T.*t3.*t9.*t11.*y+mu.*t3.*t11.*t12.*y.*3.0,-T.*t3.*t9.*t11.*z+mu.*t3.*t11.*t12.*z.*3.0,0.0,0.0,0.0,-T.*t4.*t9.*t13.*x+mu.*t4.*t12.*t13.*x.*3.0,-T.*t4.*t9.*t13.*y+mu.*t4.*t12.*t13.*y.*3.0,t15-mu.*t9-T.*t4.*t9.*t13.*z+mu.*t4.*t12.*t13.*z.*3.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[6,6]);