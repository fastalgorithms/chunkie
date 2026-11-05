function [s0,sn] = latticecoefs(n,d,kappa,l)
%CHNK.LAP2DQUAS.LATTICECOEFS precompute lattice sum coefficients for the
% quasiperiodic Laplace Green's function.
%
% Computes coefficients S_n such that the far-field periodic images of the
% quasiperiodic Laplace Green's function (images with |m| >= l) are
% represented as
%
%   G_far(x,y) = s0 + sum_{n=1}^{N} sn(n) * r^n * cos(n * theta)
%
% where r = |x-y| and theta = arg(x-y). The constant s0 is computed by
% evaluating the dual sum at a reference point and correcting with the
% explicitly summed images.
%
% Syntax: [s0,sn] = chnk.lap2dquas.latticecoefs(n,d,kappa)
%         [s0,sn] = chnk.lap2dquas.latticecoefs(n,d,kappa,l)
%
% Input:
%   n     - (N,1) vector of positive orders, e.g. 1:N
%   d     - period (scalar)
%   kappa - (nkappa,1) array of quasiperiodic phase parameters
%   l     - (optional, default 3) index of the first periodic copy
%               included in the lattice sum, i.e. copies with |m| >= l
%               are summed
%
% Output:
%   s0 - (nkappa,1) constant term of the local expansion
%   sn - (nkappa, N) polynomial coefficients for orders 1...N
%
% see also CHNK.LAP2DQUAS.GREEN, CHNK.LAP2DQUAS.KERN,
%          CHNK.LAP2DQUAS.QUASI_DUAL_SUM, CHNK.LAP2DQUAS.POLYLOG

if nargin < 4, l = 3; end

% number of copies summed explicitly on each side
lexp = l-1;

n = n(:).';
nkappa = length(kappa);

% sum_{|j| >= l} j^-n exp(i kappa d j)
%   = Li_n(e^{i kappa d}) + (-1)^n Li_n(e^{-i kappa d})
%       - sum_{0 < |j| < l} j^-n exp(i kappa d j)
sn = zeros(nkappa,length(n));
for i = 1:nkappa
    z = exp(1i*kappa(i)*d);
    sn(i,:) = chnk.lap2dquas.polylog(n,z) + (-1).^n.*chnk.lap2dquas.polylog(n,1/z);
    for j = [-lexp:-1, 1:lexp]
        sn(i,:) = sn(i,:) - j.^-n.*z^j;
    end
    sn(i,:) = sn(i,:)./n./d.^n;
end
sn = sn / 2 / pi ;
s0 = chnk.lap2dquas.quasi_dual_sum(0,d/2,0,kappa,d) ...
    - chnk.lap2dquas.green([0;0],[0;d/2],kappa,d,0,sn,lexp,1);
s0 = s0(:);

end
