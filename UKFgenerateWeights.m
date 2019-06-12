function [W0_m,W0_c,Wi_m,Wi_c] = UKFgenerateWeights(n,Beta,alpha,k)

lambda = alpha^2*(n+k)-n;
W0_m = lambda/(n+lambda);
W0_c = lambda/(n+lambda) + (1-alpha^2+Beta);
Wi_m = 1/(2*(n+lambda));
Wi_c = 1/(2*(n+lambda));