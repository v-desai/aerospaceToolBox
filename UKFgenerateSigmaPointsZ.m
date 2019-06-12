function [sigmaPointsZ] = UKFgenerateSigmaPointsZ(sigmaPoints,measurementType)

N = size(sigmaPoints,1);
estRA = zeros(1,2*N+1);
estDEC = zeros(1,2*N+1);


if measurementType == 'RADEC'
    m = 2;
    sigmaPointsZ = zeros(m,2*N+1);
    for ii = 1:(2*N+1)
        if norm(sigmaPoints(1,ii)) == 0 && norm(sigmaPoints(2,ii)) == 0
            estRA = atan2(sigmaPoints(5,ii),sigmaPoints(4,ii));
        else
            estRA = atan2(sigmaPoints(2,ii),sigmaPoints(1,ii));
            estDEC = asin(sigmaPoints(3,ii)/norm(sigmaPoints(1:3,ii)));
            sigmaPointsZ(1:m,ii) = [estRA;estDEC];
        end
    end
elseif measurementType == 'RANGE'
    rho = 0;
end