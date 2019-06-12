function [xbar] = UKFmeanStateEstimate(W0_m,Wi_m,sigmaPoints)

% Uses weights to get posterior mean and covariance
n = size(sigmaPoints,1);
xbar = zeros(n,1);
N = size(sigmaPoints,2);
for ii = 1:N
    if ii == 1
        xbar = xbar + W0_m*sigmaPoints(:,ii);
    else
        xbar = xbar + Wi_m*sigmaPoints(:,ii);
    end
end