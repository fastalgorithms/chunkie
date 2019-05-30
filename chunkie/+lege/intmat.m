function [aint,x,w] = intmat(n)
%INTMAT returns the spectral differentiation matrix on n Gaussian nodes
% a transcription of part of the Rokhlin routine legeinmt
%
% input: 
%   n - the number of Gaussian nodes 
% 

[x,w,u,v] = lege.exps(n);

coeffs = eye(n);
polints = lege.intpol(coeffs);
aint = v*polints(1:end-1,:)*u;

end
