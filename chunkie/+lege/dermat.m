function mat = dermat(k,u,v)
%LEGE.DERMAT returns the spectral differentiation matrix on Legendre nodes
% of order k
% 
% input:
%   k - integer, number of Legendre nodes
%
% optional inputs:
%   u - k x k matrix mapping values at legendre nodes to legendre series
%      coefficients
%   v - k x k matrix mapping legendre series coefficients to values at 
%      legendre nodes
%
% output:
%   mat - spectral differentiation matrix on Legendre nodes
%
%see also LEGE.EXPS, LEGE.DERPOL

if nargin < 3
    [~,~,u,v] = lege.exps(k);
end

mat = v(:,1:end-1)*lege.derpol(u);
