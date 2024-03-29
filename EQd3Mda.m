function d3Mda = EQd3Mda(a,mu,t0,tf)
%EQD3MDA
%    D3MDA = EQD3MDA(A,MU,T0,TF)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    06-Mar-2019 12:47:39

t2 = 1.0./a.^5;
t3 = mu.*t2;
t4 = t0-tf;
d3Mda = 1.0./a.^12.*mu.^2.*1.0./t3.^(3.0./2.0).*t4.*(-7.5e1./8.0)+1.0./a.^7.*mu.*1.0./sqrt(t3).*t4.*(4.5e1./2.0);
