function [sigmaPoints] = UKFgenerateSigmaPoints(xbar,Pbar,alpha,k)

n = length(xbar);
lambda = alpha^2*(n+k)-n;
L = chol(Pbar,'lower');

sigmaPoints(:,1) = xbar;
sigmaPoints(:,2:n+1) = xbar + sqrt(n+lambda)*L;
sigmaPoints(:,n+2:2*n+1) = xbar - sqrt(n+lambda)*L;

