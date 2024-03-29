function a_TwoBody = EQTwoBodyEOM(mu,x,y,z)
%EQTWOBODYEOM
%    A_TWOBODY = EQTWOBODYEOM(MU,X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    25-Feb-2019 13:08:53

t2 = x.^2;
t3 = y.^2;
t4 = z.^2;
t5 = t2+t3+t4;
t6 = 1.0./t5.^(3.0./2.0);
a_TwoBody = [-mu.*t6.*x;-mu.*t6.*y;-mu.*t6.*z];
