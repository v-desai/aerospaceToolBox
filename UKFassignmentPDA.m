function [p0,pj,zTilde,zTildej,M] = UKFassignmentPDA(measurements,zBar,PBarzz,lambdac,Pd,Pg,m,kk)

M = size(measurements,2);
dMD = zeros(M,1);
for jj = 1:M
    dMD(jj,1) = sqrt( ( measurements(:,jj,kk) - zBar )'  / PBarzz * (measurements(:,jj,kk) - zBar) );
end
b = lambdac*(1-Pd*Pg)*(2*pi)^(m/2) / Pd * sqrt(det(PBarzz));
alphaj = zeros(jj,1);
for jj = 1:M
    alphaj(jj,1) = exp(-dMD(jj,1)^2/2);
end
sumalphaj = sum(alphaj);
pj = zeros(M,1);
for jj = 0:M
    if jj == 0
        p0 = b/(b + sumalphaj);
    elseif (1<=jj) && (jj<=M)
        pj(jj) = alphaj(jj)/(b + sumalphaj);
    end
end
zTildej = zeros(m,M);
zTilde = zeros(m,1);
for jj = 1:M
    zTildej(:,jj) = measurements(:,jj,kk) - zBar;
end
for jj = 1:M
    zTilde = zTilde + pj(jj)*zTildej(:,jj);
end