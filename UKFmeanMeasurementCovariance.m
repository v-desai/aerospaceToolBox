function [Pbarzz] = UKFmeanMeasurementCovariance(W0_c,Wi_c,zbar,estZ)
% Measurement Covariance
Pbarzz = zeros(2);

for ii = 1:size(estZ,2)
    if ii == 1
        Pbarzz = Pbarzz + W0_c*(estZ(:,ii)-zbar)*(estZ(:,ii)-zbar)';
    else
        Pbarzz = Pbarzz + Wi_c*(estZ(:,ii)-zbar)*(estZ(:,ii)-zbar)';
    end
end