function S = latticecoefs(n,zk,d,kappa,alpha,a,M,l)
%CHNK.HELM2DQUAS.LATTICECOEFS precompute lattice sum coefficients for the
% quasiperiodic Helmholtz Green's function.
%
% Computes the local expansion coefficients S_n such that the far-field
% (periodic images beyond radius l) of the quasiperiodic Green's function
% is represented as
%
%   G_far(x,y) := sum_{|m|>l} i/4 H_0^{(1)}(zk |x - m d e_1 - y|)
%                   exp(i kappa d m)
%              = sum_{n=0}^{N} S_n J_n(zk r) cos(n theta)
%
% where r = |x-y| and theta = arg(x-y). The S_n are computed via
% trapezoid rule on [0,a] with M points (see Yasumoto & Yoshitomi 1999).
%
% Syntax: S = chnk.helm2dquas.latticecoefs(n,zk,d,kappa,alpha,a,M,l)
%
% Input:
%   n     - (N+1,1) vector of orders, typically 0:N
%   zk    - complex number, Helmholtz wavenumber
%   d     - period
%   kappa - (nkappa,1) array of quasiperiodic phase parameters
%   alpha - (nkappa,1) array exp(1i*kappa*d), Bloch factors
%   a     - upper limit of trapezoid quadrature interval (default 15)
%   M     - number of trapezoid quadrature nodes (default 1e4)
%   l     - number of explicit periodic copies on each side to exclude
%               (default 2)
%
% Output:
%   S - (nkappa, N+1) lattice sum coefficients; S(i,n+1) is the n-th
%       order coefficient for kappa(i)
%
% see also CHNK.HELM2DQUAS.GREEN, CHNK.HELM2DQUAS.KERN

if nargin < 8, l = 2; end
taus = linspace(0,a,M)*(1-1i);

nkappa = length(kappa);
S = zeros(nkappa,length(n));

sq_taus = sqrt(1-taus.^2);

Gp = exp(n*log(taus - 1i*sq_taus) +l*1i*zk*d*sq_taus);
Gm = exp(n*log(-taus - 1i*sq_taus) +l*1i*zk*d*sq_taus);


for i = 1:nkappa

Fp = 1./(sq_taus.*(1-exp(1i*zk*d.*(sq_taus-kappa(i)*d./zk/d))));
Fm = 1./(sq_taus.*(1-exp(1i*zk*d.*(sq_taus+kappa(i)*d./zk/d))));

S_p = sum((Gp+Gm).*Fp,2)-(Gp(:,1)+Gm(:,1)).*Fp(:,1)/2;
S_m = sum((Gp+Gm).*Fm,2)-(Gp(:,1)+Gm(:,1)).*Fm(:,1)/2;

S(i,:) = -1i*exp(1i*pi/4)*sqrt(2)/pi*((-1).^n.*alpha(i)^(-l).*S_p+alpha(i)^(l)*S_m)*a/(M-1);
end

end