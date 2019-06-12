function [componentsw,componentsx,componentsP] = GMSplitter(xBar,PBar,splitNumber,w)

[V,D] = eig(PBar);
N = size(xBar,1);
componentsw = zeros(splitNumber,1);
componentsx = zeros(N,splitNumber);
componentsP = zeros(N,N,splitNumber);

if splitNumber == 3
    alpha_ii = [0.2252246249 ; 0.5495507502 ; 0.2252246249];
    m_ii = [-1.0575154615 ; 0 ; 1.0575154615];
    sigma_ii = [0.6715662887 ; 0.6715662887 ; 0.6715662887];
elseif splitNumber == 5
    alpha_ii = [0.0763216491  ; 0.2474417860  ; 0.3524731300 ; 0.2474417860 ; 0.0763216491];
    m_ii     = [-1.6899729111 ; -0.8009283834 ; 0            ; 0.8009283834 ; 1.6899729111];
    sigma_ii = [0.4422555386  ; 0.4422555386  ; 0.4422555386 ; 0.4422555386 ; 0.4422555386];
end

index = min(find(diag(D) == max(diag(D))));
lambda_k = D(index,index);
v_k = V(:,index);
for ii = 1:splitNumber
    componentsw(ii) = alpha_ii(ii)*w;
    componentsx(:,ii) = xBar + sqrt(lambda_k)*m_ii(ii)*v_k;
    
    D_i = diag(D);
    D_i(index) = sigma_ii(ii)^2*D_i(index);
    D_ii = diag(D_i);
    
    componentsP(:,:,ii) = V*D_ii*V';
end