function dRdOE = EQdRdOE_Analytical_J2full(J2,M,R,a,e,i,mu,omega)
%EQDRDOE_ANALYTICAL_J2FULL
%    DRDOE = EQDRDOE_ANALYTICAL_J2FULL(J2,M,R,A,E,I,MU,OMEGA)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    20-Nov-2018 22:45:45

t2 = sin(i);
t3 = e.^2;
t4 = omega.*2.0;
t5 = t3.^2;
t6 = t2.^2;
t7 = M.*2.0;
t8 = R.^2;
t9 = t6.*(3.0./4.0);
t10 = t9-1.0./2.0;
t11 = cos(t7);
t12 = cos(M);
t13 = -t3+1.0;
t14 = M+t4;
t15 = cos(t14);
t16 = t4+t7;
t17 = cos(t16);
t18 = M.*4.0;
t19 = t4+t18;
t20 = cos(t19);
t21 = M.*3.0;
t22 = t4+t21;
t23 = cos(t22);
t24 = 1.0./a.^3;
t25 = e.*(7.0./2.0);
t50 = e.*t3.*(1.23e2./1.6e1);
t26 = t25-t50;
t27 = J2.*t23.*t26;
t28 = t5.*(1.3e1./1.6e1);
t51 = t3.*(5.0./2.0);
t29 = t28-t51+1.0;
t30 = J2.*t17.*t29;
t31 = t3.*(1.7e1./2.0);
t32 = t5.*(1.15e2./6.0);
t33 = t31-t32;
t34 = J2.*t20.*t33;
t35 = e./2.0;
t52 = (e.*t3)./1.6e1;
t36 = t35-t52;
t37 = t27+t30+t34-J2.*t15.*t36;
t38 = cos(i);
t39 = 1.0./t13.^(3.0./2.0);
t40 = J2.*t39;
t41 = t3.*(9.0./2.0);
t42 = t5.*(7.0./2.0);
t43 = t41+t42;
t44 = J2.*t11.*t43;
t45 = e.*3.0;
t46 = e.*t3.*(2.7e1./8.0);
t47 = t45+t46;
t48 = J2.*t12.*t47;
t49 = t40+t44+t48;
t53 = sin(t22);
t54 = sin(t16);
t55 = J2.*t29.*t54.*2.0;
t56 = sin(t19);
t57 = sin(t14);
dRdOE = [1.0./a.^4.*mu.*t8.*(t6.*t37.*(3.0./4.0)-t10.*t49).*-3.0,mu.*t8.*t24.*(t6.*(J2.*t15.*(t3.*(3.0./1.6e1)-1.0./2.0)-J2.*t23.*(t3.*(3.69e2./1.6e1)-7.0./2.0)-J2.*t17.*(e.*5.0-e.*t3.*(1.3e1./4.0))+J2.*t20.*(e.*1.7e1-e.*t3.*(2.3e2./3.0))).*(3.0./4.0)-t10.*(J2.*t12.*(t3.*(8.1e1./8.0)+3.0)+J2.*t11.*(e.*9.0+e.*t3.*1.4e1)+J2.*e.*1.0./t13.^(5.0./2.0).*3.0)),mu.*t8.*t24.*(t2.*t37.*t38.*(3.0./2.0)-t2.*t38.*t49.*(3.0./2.0)),mu.*t6.*t8.*t24.*(t55+J2.*t26.*t53.*2.0+J2.*t33.*t56.*2.0-J2.*t36.*t57.*2.0).*(-3.0./4.0),0.0,mu.*t8.*t24.*(t10.*(J2.*t43.*sin(t7).*2.0+J2.*t47.*sin(M))-t6.*(t55+J2.*t26.*t53.*3.0+J2.*t33.*t56.*4.0-J2.*t36.*t57).*(3.0./4.0))];
