function output = stabilityCheck(X)
% X is vector of real eigenvalue components
check1 = 0;
check2 = 0;
check3 = 0;

for ii = 1:length(X)
    if X(ii) < 1e-14 % eigenvalue is negative
        check1(ii) = 1;
    end
end
for ii = 1:length(X)
    if X(ii) > 1e-14 % eigenvalue is positive
        check2 = 1;
    end
end
for ii = 1:length(X)
    if X(ii) > 1e-14 && X(ii) < 1e-14 % eigenvalue is zero
        check3 = 1;
    end
end

if check1(:) == 1
    output = 'stable';
elseif check2 == 1 
    output = 'unstable';
elseif check3 == 1
    output = 'linearly stable';
end

