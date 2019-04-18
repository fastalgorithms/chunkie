function [mat,x,w,u,v] = matrin(n,xs)
%MATRIN form the matrix for interpolating from n Legendre nodes to
% the points xs
%
% optionally returns the roots, weights, and matrices for switching
% between values and coefficients (u : values -> coeffs, 
% v : coeffs -> values)

[x,w,u,v] = lege.exps(n);

xsvals = lege.pols(xs,n-1);

mat = (xsvals).'*u;

end