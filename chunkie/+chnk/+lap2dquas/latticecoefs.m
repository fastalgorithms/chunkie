function [s0,sn] = latticecoefs(n,d,kappa,l)
%CHNK.LAP2DQUAS.LATTICECOEFS precompute lattice sums for the quasi-periodic
% Laplace problem
%
% quasi periodic parameters:
%   d - period
%   kappa - quasiperiodic parameters
%   l - number of periodic copies on each side to exclude

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