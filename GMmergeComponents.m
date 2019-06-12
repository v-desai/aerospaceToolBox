function [weights,xHatNew,PNew] = GMmergeComponents(pruningThreshhold,mergingThreshhold,maxComponents,wHat,xHatGM,PGM)

% Orbital Debris - Lecture 5 Slide 12 & 13

N = size(xHatGM,1);
% Step 1
I = find(wHat >= pruningThreshhold);
current = 0;
% Step 2
while ~isempty(I)
    % Step 2.1
    current = current + 1;
    % Step 2.2
    jj = find(wHat(I) == max(wHat(I)));
    jj = I(min(jj));
    % Step 2.3
    dMDsquared = zeros(1,length(I));
    for ii = I
         dMDsquared(1,ii) = (xHatGM(1:6,ii) - xHatGM(1:6,jj))' / (PGM(1:6,1:6,ii)) * (xHatGM(1:6,ii)-xHatGM(1:6,jj)) ;
    end
    L = intersect(I,find(dMDsquared <= mergingThreshhold^2));
    % Step 2.4
    wHatMerge = sum(wHat(L));
    xHatMerge = zeros(N,1);
    for ii = L
        xHatMerge = xHatMerge + wHat(ii)/wHatMerge*xHatGM(:,ii);
    end
    PMerge = zeros(N,N);
    for ii = L
        PMerge = PMerge + wHat(ii)/wHatMerge*( PGM(:,:,ii) + xHatGM(:,ii)*xHatGM(:,ii)' );
    end
    PMerge = PMerge - xHatMerge*xHatMerge';
    
    weightsNew(:,current) = wHatMerge;
    xHatNew(:,current) = xHatMerge;
    PNew(:,:,current) = PMerge;
    % Step 2.5
    I = setdiff(I,L);
end
% Step 3
if current > maxComponents
    weightsNew = weightsNew(:,1:maxComponents);
    xHatNew = xHatNew(:,1:maxComponents);
    PNew = PNew(:,:,1:maxComponents);thi    
end
% Step 4
weights = zeros(1,current);
for ii = 1:current
    weights(ii) = weightsNew(ii)/sum(weightsNew);
end
