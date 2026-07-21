function [s0,sn] = latticecoefs(n,d,kappa,l)
%CHNK.LAP2DQUAS.LATTICECOEFS precompute lattice sum coefficients for the
% quasiperiodic Laplace Green's function.
%
% Computes coefficients S_n such that the far-field periodic images of the
% quasiperiodic Laplace Green's function (beyond radius l) are represented as
%
%   G_far(x,y) = s0 + sum_{n=1}^{N} sn(n) * r^n * cos(n * theta)
%
% where r = |x-y| and theta = arg(x-y). The constant s0 is computed by
% evaluating the dual sum at a reference point and correcting with the
% known free-space Green's function.
%
% Syntax: [s0,sn] = chnk.lap2dquas.latticecoefs(n,d,kappa)
%         [s0,sn] = chnk.lap2dquas.latticecoefs(n,d,kappa,l)
%
% Input:
%   n     - (N,1) vector of positive orders, e.g. 1:N
%   d     - period (scalar)
%   kappa - (nkappa,1) array of quasiperiodic phase parameters
%   l     - (optional, default 2) number of explicit periodic copies on
%               each side to exclude from the lattice sum
%
% Output:
%   s0 - (nkappa,1) constant term of the local expansion
%   sn - (nkappa, N) polynomial coefficients for orders 1..N
%
% see also CHNK.LAP2DQUAS.GREEN, CHNK.LAP2DQUAS.KERN,
%          CHNK.LAP2DQUAS.QUASI_DUAL_SUM

if nargin < 4, l = 2; end

nkappa = length(kappa);


sn = 0;
for j = [-l:-1, 1:l]
    sn = sn + j.^-n.*exp(1i*kappa(:)*d * j);
end

s0 = zeros(nkappa,1);
for i = 1:nkappa
    sn(i,:) = -sn(i,:) + polylog(n,exp(1i*kappa(i)*d)) + (-1).^-n.*polylog(n,exp(-1i*kappa(i)*d));
    sn(i,:) = sn(i,:)./n./d.^n;
end
sn = sn / 2 / pi ;
s0 = chnk.lap2dquas.quasi_dual_sum(0,d/2,0,kappa,d)- chnk.lap2dquas.green([0;0],[0;d/2],kappa,d,0,sn,l,1);




end