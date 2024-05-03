function [aint,u,v] = intmat(n,u,v)
%INTMAT returns the spectral integration matrix on n Gaussian nodes
% a transcription of part of the Rokhlin routine legeinmt
%
% input: 
%   n - the number of Gaussian nodes 
%
% optional inputs:
%   u - k x k matrix mapping values at legendre nodes to legendre series
%      coefficients
%   v - k x k matrix mapping legendre series coefficients to values at 
%      legendre nodes
% 
% output: 
%   

if nargin < 3
    [~,~,u,v] = lege.exps(n);
end

tmp = lege.intpol(u,'original');
aint = v*tmp(1:end-1,:);

end
