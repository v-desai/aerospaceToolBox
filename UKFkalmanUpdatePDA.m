function [xHat,P] = UKFkalmanUpdatePDA(pj,p0,PBarxz,PBarzz,xBar,PBar,zTilde,zTildej,Q,n,d,N,M)

K = PBarxz/PBarzz;
xHat = xBar + K*zTilde;
xHat(n+1:N,1) = zeros(d,1);

sumpjzj = 0;
for jj = 1:M
    sumpjzj = sumpjzj + pj(jj)*zTildej(:,jj)*zTildej(:,jj)';
end

P = p0*PBar + (1-p0)*(PBar - K*PBarzz*K') + K*(sumpjzj-zTilde*zTilde')*K';
Pnew = zeros(N); 
Pnew(1:n,1:n) = P(1:n,1:n); 
Pnew(n+1:N,n+1:N) = Q;
P = Pnew;