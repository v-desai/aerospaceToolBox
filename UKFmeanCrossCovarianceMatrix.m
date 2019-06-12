function [Pbarxz] = UKFmeanCrossCovarianceMatrix(W0_c,Wi_c,sigmaPoints,xbar,estZ,zbar)

% Cross-covariance Matrix
n = size(sigmaPoints,1);
m = size(estZ,1);
Pbarxz = zeros(n,m);

for ii = 1:(2*n+1)
    if ii == 1
        Pbarxz = Pbarxz + W0_c*(sigmaPoints(:,ii)-xbar)*(estZ(1:m,ii)-zbar)';
    else
        Pbarxz = Pbarxz + Wi_c*(sigmaPoints(:,ii)-xbar)*(estZ(1:m,ii)-zbar)';
    end
end