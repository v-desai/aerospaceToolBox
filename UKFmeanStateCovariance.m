function [Pbar] = UKFmeanStateCovariance(W0_c,Wi_c,xbar,sigmaPoints)

% Uses weights to get posterior mean and covariance
n = size(sigmaPoints,1);
Pbar = zeros(n);

for ii = 1:(2*n+1)
    if ii == 1
        Pbar = Pbar + W0_c*(sigmaPoints(:,ii)-xbar)*(sigmaPoints(:,ii)-xbar)';
    else
        Pbar = Pbar + Wi_c*(sigmaPoints(:,ii)-xbar)*(sigmaPoints(:,ii)-xbar)';
    end
end
