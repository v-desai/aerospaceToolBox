function [assignment] = measurementAssignment(J,observations,zbarGM,PbarzzGM,Pg)

costTable = zeros(J,size(observations,2));

for ii = 1:J
    for jj = 1:size(observations,2)
        dMD = ( observations(:,jj) - zbarGM(:,ii) )'  / PbarzzGM(:,:,ii) * (observations(:,jj) - zbarGM(:,ii)) ;
        costTable(ii,jj) = dMD;
    end
end

G = chi2inv(Pg,2) ^ 0.5;
birthCost = Inf(size(observations,2));
for ii = 1:size(birthCost,2)
    birthCost(ii,ii) = G^2;
end

costMat = [costTable ; birthCost];
[assignment,~] = munkres(costMat);