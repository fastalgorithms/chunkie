function [mat,x,w,u,v] = matrin(n,xs,varargin)
%LEGE.MATRIN form the matrix for interpolating from n Legendre nodes to
% the points xs
%
% optionally returns the roots, weights, and matrices for switching
% between values and coefficients (u : values -> coeffs, 
% v : coeffs -> values), which is the output of LEGE.EXPS

if nargin > 2
    u = varargin{1};
end

if (nargin < 3 || nargout > 1)

    [x,w,u,v] = lege.exps(n);
end

xsvals = lege.pols(xs,n-1);

mat = (xsvals).'*u;

end