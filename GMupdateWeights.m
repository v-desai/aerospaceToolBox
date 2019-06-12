function [wHat] = GMupdateWeights(w,zGM,zBarGM,PBarzzGM)

% Orbital Debris: Lecture 5 - Slide 7

J = size(zGM,2);
wTilde = zeros(1,J);

for ii = 1:J
    wTilde(ii) = w(ii)*mvnpdf(zGM(:,ii),zBarGM(:,ii),PBarzzGM(:,:,ii));
end

wHat = wTilde./sum(wTilde);