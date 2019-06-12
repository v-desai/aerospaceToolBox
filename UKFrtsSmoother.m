% Unscented RTS Fixed-Interval Smoother
function [xS,PS,xpost,time] = UKFrtsSmoother(filename)
% ---------------------------------------------------------------------
% Description:
%
%  Takes output of UKF and smooths the data
%
% Inputs:
%
%  filename - name of file containing below variables
%
%  time  - vector of time of measurement
%  xpre  - apriori state means
%  xpost - aposteriori state means
%  Ppre  - apriori covariance
%  Ppost - aposteriori covariance
%  SPpre - sigma points without process noise applied
%
% Outputs:
%
%  xS - Smoothed state matrix
%  PS - Smoothed covariance matrix
%  time - time list

% Assumptions/References:
%
%  List papers Saarka
%
% Dependencies:
%
%  UKFgenerateWeights
%  UKFgenerateSigmaPoints
%
% Modification History:
%
%  19march19     Vivek Desai      original version
%
% ---------------------------------------------------------------------
% Copyright ASTRIA, 2019
% ---------------------------------------------------------------------

%% ========== Load Data ==========
jsonCheck = filename(end-3:end);
if jsonCheck == 'json'
    jsonData = jsondecode(fileread(filename));
    for ii = 1:length(jsonData.Estimation)
        xpre(:,ii) = jsonData.Estimation(ii).xBar;
        xpost(:,ii) = jsonData.Estimation(ii).EstimatedState;
        Ppre(:,:,ii) = jsonData.Estimation(ii).PBar;
        Ppost(:,:,ii) = jsonData.Estimation(ii).EstimatedCovariance;
    end
    T = size(xpost,2);
    time = 1:length(jsonData.Estimation);       
    N = size(xpost,1);
    Wi_c = 1/(2*N);
    Beta = 2; alpha = 1; k = 3-N;
else
    load(filename)    
    xpre = store_xBar;
    xpost = store_xHat;
    Ppre = store_PBar;
    Ppost = store_P;
    T = size(xpost,2);
    time = 1:T;
    N = size(xpost,1);
    Beta = 2; alpha = 1; k = 3-N;
    [~,~,~,Wi_c] = UKFgenerateWeights(N,Beta,alpha,k);
end

%% ========== Run Smoother ==========
for kk = T:-1:1
    
    if kk == T
        xS(:,kk) = xpost(:,kk);
        PS(:,:,kk) = Ppost(:,:,kk);
    else
        [SPpre] = UKFgenerateSigmaPoints(xpre(:,kk+1),Ppre(:,:,kk+1),alpha,k);
        [SPpost] = UKFgenerateSigmaPoints(xpost(:,kk),Ppost(:,:,kk),alpha,k);
        C = zeros(N);
        for ii = 1:2*N+1
            C = C + Wi_c*(SPpost(:,ii)-xpost(:,kk))*(SPpre(:,ii)-xpre(:,kk+1))';
        end
        A = C/Ppre(:,:,kk+1);
        xS(:,kk) = xpost(:,kk) + A*(xS(:,kk+1)-xpre(:,kk+1));
        PS(:,:,kk) = Ppost(:,:,kk) + A*(PS(:,:,kk+1)-Ppre(:,:,kk+1))*A';
    end
    
end

stateDiff = xS - xpost;
covDiff = PS - Ppost;
for ii = 1:size(stateDiff,1)
    sigma = sqrt(covDiff(1,1));
    constistencyRatio = stateDiff./sigma;
end

a = 0;


