% rvPropTEST

%% TEST rvProp function

%% ---------- Two-Body + J2 + J3 -------------------------------------
% Initial conditions and solution from Dr. Jones - Orbital Debis HW 1

x0 = [-2012.151 ; -381.450 ; 6316.615 ; 5.400336 ; -5.916814 ; 1.362965];
xfTrue = [ -2006.89250117544 ; 4499.62892977502 ; -6106.21106747767 ; -4.4332289059933 ; 2.98092880414086 ; 4.37202659137548 ];

t0 = 0;
tf = 86400;
dt = 30;

runTwoBody = 1;
runJ2 = 1;
runJ3 = 1;
runDrag = 0;
runSTM = 0;
runUKF = 0;
runAEGIS = 0;
input = [runTwoBody,runJ2,runJ3,runDrag,runSTM,runUKF,runAEGIS];
options = odeset('reltol',1e-13,'abstol',1e-13);

[~,RV] = ode45(@rvProp,t0:dt:tf,x0,options,input);

xf = RV(end,:)';

error = xf - xfTrue;

fprintf('%.11f\n',error(1:3))
fprintf('%.13f\n',error(4:6))

fprintf('\n')
%% ---------- Two-Body + J2 + J3 + Drag -------------------------------------
% Initial conditions and solution from Dr. Jones - Orbital Debis HW 1
% x0 = [-2012.151 ; -381.450 ; 6316.615 ; 5.400336 ; -5.916814 ; 1.362965];
% xfTrue = [-2013.38695085811 ; 4503.89575032041 ; -6099.5374817889 ; -4.430906534889688 ; 2.97557481433479 ; 4.37949674436265];

% 2nd Set - Orbital Debris HW 2 & 3
x0 = [-2011.990 ; -382.065 ; 6316.376 ; 5.419783 ; -5.945319 ; 1.37398];
xfTrue = [4105.111 ; -1610.374 ; -6985.989 ; -3.241787 ; 4.690008 ; -3.257038];

t0 = 0;
tf = 86400;
dt = 30;

runTwoBody = 1;
runJ2 = 1;
runJ3 = 1;
runDrag = 1;
runSTM = 0;
runUKF = 0;
runAEGIS = 0;
input = [runTwoBody,runJ2,runJ3,runDrag,runSTM,runUKF,runAEGIS];
options = odeset('reltol',1e-13,'abstol',1e-13);

[~,RV] = ode45(@rvProp,t0:dt:tf,x0,options,input);

xf = RV(end,:)';

error = xf - xfTrue;

% fprintf('%.11f\n',error(1:3))
% fprintf('%.13f\n',error(4:6))
fprintf('%.3f\n',error(1:3))
fprintf('%.6f\n',error(4:6))

%% -----------------------------------------------------------------------

x0 = [-2011.990 ; -382.065 ; 6316.376 ; 5.419783 ; -5.945319 ; 1.37398];
P0 = [eye(3),zeros(3) ; zeros(3),10^-6*eye(3)];

t0 = 0;
tf = 86400;
dt = 30;

runTwoBody = 1;
runJ2 = 1;
runJ3 = 1;
runDrag = 1;
runSTM = 0;
runUKF = 0;
runAEGIS = 0;
input = [runTwoBody,runJ2,runJ3,runDrag,runSTM,runUKF,runAEGIS];
options = odeset('reltol',1e-13,'abstol',1e-13);

[componentsw,componentsxBar,componentsPBar] = GMSplitter(x0,P0,5);

for ii = 1:5
[T,RVset] = ode45(@rvProp,t0:dt:tf,componentsxBar(:,ii),options,input);
xf(:,ii) = RVset(end,:)';
end

disp(xf)