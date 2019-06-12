function [xHat,P] = UKFkalmanUpdate(PBarxz,PBarzz,xBar,PBar,z,zBar,Q,n,d,N)

K = PBarxz/PBarzz;
xHat = xBar + K*(z - zBar);
xHat(n+1:N,1) = zeros(d,1);
P = PBar - K*PBarzz*K';
Pnew = zeros(N);
Pnew(1:n,1:n) = P(1:n,1:n); 
Pnew(n+1:N,n+1:N) = Q;
P = Pnew;