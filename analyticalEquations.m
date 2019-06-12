%% ANALYTICAL EQUATIONS
clear all; clc;
filePath = [userpath '\aerospaceToolboxVD'];
addpath(filePath);
cd(filePath)
equationSet = 15;
% case 1 : acceleration (Two Body + J2)
% case 2 : J2 potential function given orbital elements
% case 3 : Averaged J2 potential given orbital elements
% case 4 : dRd0E equation given J2 potential
% case 5 : dRdOE equation given average J2 potential
% case 6 : Partials of position and velocity w.r.t. orbital elements
% case 7 : Simplified Lagrange Bracket
% case 8 : acceleration (CW) *with a radial thrust
% case 9 : CW STM *with a radial thrust
% case 10: Computing A for CRTBP (acceleration Three Body Linear & NL)
% case 11: acceleration (J3)
% case 12: acceleration (Two Body)
% case 13: acceleration (J2)
% case 14: acceleration (Drag)
% case 15: A matrix
% case 16: STT 2nd order
% case 17: STT 3nd order
% case 18: STT nth order

switch equationSet
    
    case 1 % acceleration (Two Body + J2)
        syms x y z mu r R J2
        r = (x^2 + y^2 + z^2)^.5;
        U = mu/r*(1 - R^2/r^2*J2*(3*z^2/(2*r^2)-1/2));
        ax_J2 = diff(U,x);
        ay_J2 = diff(U,y);
        az_J2 = diff(U,z);
        a_J2(1:3,1) = [ax_J2;ay_J2;az_J2];
        matlabFunction(a_J2,'File','EQJ2EOM');
        
    case 2 % J2 potential function given orbital elements
        syms a e i omega M R mu J2
        t1 = a ^ 2;
        t5 = sin(i);
        t6 = t5 ^ 2;
        t8 = (e ^ 2);
        t9 = t8 * e;
        t13 = 2 * omega;
        t15 = cos((t13 + M));
        t18 = t8 ^ 2;
        t23 = cos((2 * omega + 2 * M));
        t31 = cos((t13 + 3 * M));
        t39 = cos((t13 + 4 * M));
        t51 = cos((2 * M));
        t58 = cos(M);
        t61 = 1 - t8;
        t62 = sqrt(t61);
        t69 = R ^ 2;
        U_J2 = mu / t1 / a * (-0.3e1 / 0.4e1 * t6 * (-(-e / 0.2e1 + t9 / 0.16e2) * J2 * t15 - (0.1e1 - 0.5e1 / 0.2e1 * t8 + 0.13e2 / ...
            0.16e2 * t18) * J2 * t23 - (0.7e1 / 0.2e1 * e - 0.123e3 / 0.16e2 * t9) * J2 * t31 - (0.17e2 / 0.2e1 * t8 - 0.115e3 / 0.6e1 * ...
            t18) * J2 * t39) + (0.3e1 / 0.4e1 * t6 - 0.1e1 / 0.2e1) * (-0.2e1 * (0.9e1 / 0.4e1 * t8 + 0.7e1 / 0.4e1 * t18) * J2 * t51 - ...
            0.2e1 * (0.3e1 / 0.2e1 * e + 0.27e2 / 0.16e2 * t9) * J2 * t58 - 0.1e1 / t62 / t61 * J2)) * t69;
        matlabFunction(U_J2,'File','EQU_Analytical_J2');
        
    case 3 % Averaged J2 potential given orbital elements
        syms J2 M R a e i mu omega
        U_J2 = EQU_Analytical_J2(J2,M,R,a,e,i,mu,omega);
        U_J2avg = int(U_J2,M,0,2*pi);
        U_J2avg = simplify(U_J2avg/(2*pi));
        matlabFunction(U_J2avg,'File','EQU_Analytical_J2avg');
        
    case 4 % dRd0E equation given J2 potential
        syms a e i omega Omega M J2 R mu
        OE = [a e i omega Omega M];
        U_J2 = EQU_Analytical_J2(J2,M,R,a,e,i,mu,omega);
        disturbingPotential = U_J2;
        dRdOE(1) = diff(disturbingPotential,OE(1));
        dRdOE(2) = diff(disturbingPotential,OE(2));
        dRdOE(3) = diff(disturbingPotential,OE(3));
        dRdOE(4) = diff(disturbingPotential,OE(4));
        dRdOE(5) = diff(disturbingPotential,OE(5));
        dRdOE(6) = diff(disturbingPotential,OE(6));
        matlabFunction(dRdOE,'File','EQdRdOE_Analytical_J2full');
        
    case 5 % dRdOE equation given average J2 potential
        syms a e i omega Omega M J2 R mu
        OE = [a e i omega Omega M];
        U_J2avg = EQU_Analytical_J2avg(J2,R,a,e,i,mu);
        disturbingPotential = U_J2avg;
        dRdOE(1) = diff(disturbingPotential,OE(1));
        dRdOE(2) = diff(disturbingPotential,OE(2));
        dRdOE(3) = diff(disturbingPotential,OE(3));
        dRdOE(4) = diff(disturbingPotential,OE(4));
        dRdOE(5) = diff(disturbingPotential,OE(5));
        dRdOE(6) = diff(disturbingPotential,OE(6));
        matlabFunction(dRdOE,'File','EQdRdOE_Analytical_J2avg');
        
    case 6 % Partials of position and velocity w.r.t. orbital elements
        syms a e i omega Omega M0 nu t0 t mu E real
        assume(e>0 & e<1)
        assume(a>0)
        assume(t0>=0)
        assume(t>=0)
        assume(mu>0)
        OE = [a e i omega Omega M0]';
        
        rabs = a*(1-e^2)/(1+e*cos(nu));
        h = (mu*a*(1-e^2))^(1/2);
        p=a*(1-e^2);
        
        r(1) = rabs*(cos(Omega)*cos(omega+nu) - sin(Omega)*sin(omega+nu)*cos(i));
        r(2) = rabs*(sin(Omega)*cos(omega+nu) + cos(Omega)*sin(omega+nu)*cos(i));
        r(3) = rabs*(sin(i)*sin(omega+nu));
        v(1) = (r(1)*h*e)/(rabs*p)*sin(nu) - h/rabs*(cos(Omega)*sin(omega+nu) + sin(Omega)*cos(omega+nu)*cos(i));
        v(2) = (r(2)*h*e)/(rabs*p)*sin(nu) - h/rabs*(sin(Omega)*sin(omega+nu) - cos(Omega)*cos(omega+nu)*cos(i));
        v(3) = (r(3)*h*e)/(rabs*p)*sin(nu) + h/rabs*sin(i)*cos(omega+nu);
        RV = [r v]';
        
        nu = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
        dnudE = diff(nu,E);
        dnude = diff(nu,e);
        dEdM0 = 1/(1-e*cos(E));
        dEde = sin(E)/(1-e*cos(E));
        dEda = (-3/2)*sqrt(mu/a^5)*(t-t0)/(1-e*cos(E));
        
        dnudM0 = dnudE*dEdM0;
        dnude = dnudE*dEde+dnude;
        dnuda = dnudE*dEda;
        
        dnudc = [dnuda dnude 0 0 0 dnudM0];
        
        syms nu
        dRVdOE_M0 = jacobian(RV,OE') + diff(RV,nu)*dnudc;
        matlabFunction(dRVdOE_M0,'File','EQdRVdOE_Analytical');
        
    case 7 % Simplified Lagrange Bracket
        syms mu e a i
        Lbrack(4,1) = (1/2)*mu*sqrt(1-e^2)/(sqrt(mu/a)*a);
        Lbrack(5,1) = (1/2)*mu*cos(i)*sqrt(1-e^2)/(sqrt(mu/a)*a);
        Lbrack(6,1) = (1/2)*mu/(sqrt(mu/a)*a);
        Lbrack(4,2) = -mu*e/(sqrt(1-e^2)*sqrt(mu/a));
        Lbrack(5,2) = -mu*e*cos(i)/(sqrt(1-e^2)*sqrt(mu/a));
        Lbrack(5,3) = -a*sin(i)*sqrt(1-e^2)*sqrt(mu/a);
        Lbrack(1,4) = -Lbrack(4,1);
        Lbrack(1,5) = -Lbrack(5,1);
        Lbrack(1,6) = -Lbrack(6,1);
        Lbrack(2,4) = -Lbrack(4,2);
        Lbrack(2,5) = -Lbrack(5,2);
        Lbrack(3,5) = -Lbrack(5,3);
        matlabFunction(Lbrack,'File','EQLbrack_Analytical');
        
    case 8 % acceleration (CW) *with a radial thrust
        syms x y z vx vy vz R mu ux uy uz T
        rvec = [x;y;z];
        vvec = [vx;vy;vz];
        uvec = [ux;uy;uz];
        R = norm(rvec);
        a_TwoBody = -mu*rvec/R^3 + rvec/norm(R)*T;
        matlabFunction(a_TwoBody,'File','EQa_TwoBody');
        
    case 9 % CW STM *with a radial thrust
        syms x y z vx vy vz R mu ux uy uz T
        rvec = [x;y;z];
        vvec = [vx;vy;vz];
        a_TwoBody = EQa_TwoBody(T,mu,x,y,z);
        F_TwoBody = jacobian([vvec;a_TwoBody],[rvec;vvec]);
        matlabFunction(F_TwoBody,'File','EQF_TwoBody');
        
    case 10 % Computing A for CRTBP
        syms x y z xdot ydot zdot MU Ydot0 dY
        r1 = sqrt((x+MU)^2+y^2+z^2);
        r2 = sqrt((x-1+MU)^2+y^2+z^2);
        % Full Nonlinear EOM (p362)
        F = 1/2*(x^2+y^2) + (1-MU)/r1 + MU/r2;
        Fx = diff(F,x);
        Fy = diff(F,y);
        Fz = diff(F,z);
        Y = [x;y;z;xdot;ydot;zdot];
        Ydot = [xdot;ydot;zdot;Fx+2*ydot;Fy-2*xdot;Fz];
        A = jacobian(Ydot,Y);
        matlabFunction(A,'File','EQA_CRTBP');
        matlabFunction(Ydot,'File','EQa_CRTPB_NL');
        % Linearized EOM
        dYdot = A*dY;
        Ydot = Ydot0 + dYdot;
        matlabFunction(Ydot,'File','EQa_CRTPB_Linear');
        
    case 11 % acceleration (J3)
        syms x y z mu r R J3
        r = (x^2 + y^2 + z^2)^.5;
        U = -J3*mu/(2*r)*(R/r)^3*(5*(z/r)^3-3*(z/r));
        ax_J3 = diff(U,x);
        ay_J3 = diff(U,y);
        az_J3 = diff(U,z);
        a_J3(1:3,1) = [ax_J3;ay_J3;az_J3];
        matlabFunction(a_J3,'File','EQJ3EOM');
        
    case 12 % acceleration (Two Body)
        syms x y z mu r R J2
        r = (x^2 + y^2 + z^2)^.5;
        U = mu/r;
        ax_TwoBody = diff(U,x);
        ay_TwoBody = diff(U,y);
        az_TwoBody = diff(U,z);
        a_TwoBody(1:3,1) = [ax_TwoBody;ay_TwoBody;az_TwoBody];
        matlabFunction(a_TwoBody,'File','EQTwoBodyEOM');
        
    case 13 % acceleration (J2)
        syms x y z mu r R J2
        r = (x^2 + y^2 + z^2)^.5;
        U = -3*J2*mu/(2*r)*(R/r)^2*((z/r)^2-1/3);
        ax_J2 = diff(U,x);
        ay_J2 = diff(U,y);
        az_J2 = diff(U,z);
        a_J2(1:3,1) = [ax_J2;ay_J2;az_J2];
        matlabFunction(a_J2,'File','EQJ2EOM');
        
    case 14 % acceleration (Drag)
        syms Cd A m thetadot R x y z xdot ydot zdot rho0 h0 H h
        rhoA = rho0*exp(-(h-h0)/H);
        Va = [xdot+thetadot*y ; ydot - thetadot*x ; zdot];
        a_Drag = -1/2*Cd*A/m*rhoA*norm(Va)*Va;
        matlabFunction(a_Drag,'File','EQDragEOM');
        
    case 15 % A matrix
        runTwoBody = 1;
        runJ2      = 1;
        runJ3      = 1;
        runDrag    = 0;
        a_Total = zeros(3,1);
        if runTwoBody
            syms mu x y z vx vy vz
            a_Total = a_Total + EQTwoBodyEOM(mu,x,y,z);
        end
        if runJ2
            syms J2 R
            a_Total = a_Total + EQJ2EOM(J2,R,mu,x,y,z);
        end
        if runJ3
            syms J3
            a_Total = a_Total + EQJ3EOM(J3,R,mu,x,y,z);
        end
        if runDrag
            syms A Cd H h h0 m rho0 thetadot xdot ydot zdot
            a_Total = a_Total + EQDragEOM(A,Cd,H,h,h0,m,rho0,thetadot,x,xdot,y,ydot,zdot);
        end
        RVsym = [x ; y ; z ; vx ; vy ; vz];
        dRVsym = [vx ; vy ; vz ; a_Total];
        A = jacobian(dRVsym,RVsym);
        matlabFunction(A,'File','EQAmatrixEOM');
        
    case 16 % STT 2nd order
        syms mu a tf t0
        dMda = -3/2*sqrt(mu/a^5)*(tf-t0);
        d2Mda = diff(dMda,a);
        matlabFunction(d2Mda,'File','EQd2Mda');
        
    case 17 % STT 3nd order
        syms mu a tf t0
        dMda = -3/2*sqrt(mu/a^5)*(tf-t0);
        d2Mda = diff(dMda,a);
        d3Mda = diff(d2Mda,a);
        matlabFunction(d3Mda,'File','EQd3Mda');
        
     case 18 % STT nth order
        syms mu a tf t0 p
        dMda = -3/2*sqrt(mu/a^5)*(tf-t0);
        dMda_nth = diff(dMda,0,a);
        p = 2;
        nthdiff = ((-1)^p * (3/2)^p * mu * (tf-t0))/( a^(p+3)*(mu/a^3)^(1/2));
end
