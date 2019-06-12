function [zbar] = UKFmeanMeasurementEstimate(W0_m,Wi_m,estZ)

% Mean Measurement Estimate
m = size(estZ,1);
N = size(estZ,2);
zbar = zeros(m,1);

for ii = 1:N
    if ii == 1
        zbar = zbar + W0_m*estZ(1:m,ii);
    else
        zbar = zbar + Wi_m*estZ(1:m,ii);
    end
end