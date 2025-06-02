function S = latticecoefs(n,zk,d,kappa,alpha,a,M,l)
%CHNK.HELMQUAS.LATTICECOEFS precompute lattice sum integrals 
% i.e. local coefficients for the representation
%   G_far(x,y) := sum_|n|>l i/4 H_0^{(1)}(zk |x- n d e_1-y|) 
%                   exp(i kapppa d n)
%              = sum_n=0^inf Sn J_n(zk |x-y|) cos(n arg(x-y))
% 
% (see Yasumoto and K. Yoshitomi (1999)) 
%
% computes the nth order lattice sum using trapezoid rule on the interval 
%   [0,a] with M points
%
% quasi periodic parameters:
%   d - period
%   kappa - quasiperiodic parameters
%   alpha - exp(1i*d*kappa), chosen so that u(x+d) = alpha * u(x)

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