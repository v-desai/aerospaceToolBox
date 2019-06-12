function [xBar,P] = GMupdateMeanCovariance(xHatGM,PGM,wHat)

% Orbital Debris: Lecture 5 - Slide 8

N = size(xHatGM,1);
J = size(xHatGM,2);

xBar = zeros(N,1);
P = zeros(N,N);

for i = 1:J
    xBar = xBar + wHat(i)*xHatGM(:,i);
end
for i = 1:J
    P = P + wHat(i)*( PGM(:,:,i) + (xHatGM(:,i)-xBar) * (xHatGM(:,i)-xBar)' );
end