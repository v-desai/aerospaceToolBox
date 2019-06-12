function a_TwoBody = EQa_TwoBody(T,mu,x,y,z)
%EQA_TWOBODY
%    A_TWOBODY = EQA_TWOBODY(T,MU,X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    20-Nov-2018 22:53:25

t2 = abs(x);
t3 = abs(y);
t4 = abs(z);
t5 = t2.^2;
t6 = t3.^2;
t7 = t4.^2;
t8 = t5+t6+t7;
t9 = 1.0./sqrt(t8);
t10 = 1.0./t8.^(3.0./2.0);
a_TwoBody = [T.*t9.*x-mu.*t10.*x;T.*t9.*y-mu.*t10.*y;T.*t9.*z-mu.*t10.*z];
