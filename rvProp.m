function dRV = rvProp(t,RV,input)

% options = odeset('reltol',1e-13,'abstol',1e-13);
% input = [runTwoBody,runJ2,runJ3,runDrag,runSTM,runUKF,runAEGIS,runLaser]
% [T,RVset] = ode45(@rvProp,[time(kk-1) time(kk)],RV ,options, input);

dRV = zeros(length(RV),1);  % Pre-allocate
x = RV(1);                  % [km]
y = RV(2);                  % [km]
z = RV(3);                  % [km]
dRV(1:3,1) = RV(4:6,1);     % [km/s]

%% Edit Peturbations
% Constants
mu = 398600.4415;   % [km^3/s^2] Earth Gravitational Parameter
J2 = 0.0010826267;  % Earth J2
J3 = -0.0000025327; % Earth J3
R  = 6378.1363;     % [km] Earth Radius

% TwoBody
if isfield(input,'runTwoBody')
    dRV(4:6,1) = dRV(4:6,1) + EQTwoBodyEOM(mu,x,y,z);
end

% J2
if isfield(input,'runJ2')
    dRV(4:6,1) = dRV(4:6,1) + EQJ2EOM(J2,R,mu,x,y,z);
end

% J3
if isfield(input,'runJ3')
    dRV(4:6,1) = dRV(4:6,1) + EQJ3EOM(J3,R,mu,x,y,z);
end

% Drag
if isfield(input,'runDrag')
    h = abs(R - norm(RV(1:3)));
    [rho0,h0,H] = getDensityParams(h);
    rho0 = rho0 * 1000^3; % [kg/km^3]
    xdot = RV(4);         % [km/s]
    ydot = RV(5);         % [km/s]
    zdot = RV(6);         % [km/s]
    Cd = input.runDrag(2);
    A = input.runDrag(3);
    m = input.runDrag(4);
    thetadot = input.runDrag(5);
    a_Drag = EQDragEOM(A,Cd,H,h,h0,m,rho0,thetadot,x,xdot,y,ydot,zdot);
    dRV(4:6,1) = dRV(4:6,1) + a_Drag;
end

% State Transition Matrix
if isfield(input,'runSTM')
    STM = reshape(RV(7:42),6,6);
    Amatrix = EQAmatrixEOM(A,Cd,H,J2,J3,R,h,h0,m,mu,rho0,thetadot,x,xdot,y,ydot,z,zdot);
    STMdot = Amatrix*STM;
    dRV(7:42,1) = reshape(STMdot,36,1);
end

if isfield(input,'runUKF')
    dRV(4:6,1) = dRV(4:6,1) + RV(7:9,1); % Add process noise
    dRV(7:9,1) = 0; % Derivative of process noise and measurement noise
end

if isfield(input,'runAEGIS')
    
    for ii = 4:6:76
        
        x = RV(ii-3);                  % [km]
        y = RV(ii-2);                  % [km]
        z = RV(ii-1);                  % [km]
        dRV(ii-3:ii-1,1) = RV(ii:ii+2,1);
        dRV(ii:ii+2,1) = dRV(ii:ii+2,1) + EQTwoBodyEOM(mu,x,y,z);
        dRV(ii:ii+2,1) = dRV(ii:ii+2,1) + EQJ2EOM(J2,R,mu,x,y,z);
        dRV(ii:ii+2,1) = dRV(ii:ii+2,1) + EQJ3EOM(J3,R,mu,x,y,z);
        
        if ii == 4
            Amatrix = EQAmatrixEOM(J2,J3,R,mu,x,y,z);
            entropyDot = trace(Amatrix);
            dRV(79,1)= entropyDot;
        end
        
    end

end