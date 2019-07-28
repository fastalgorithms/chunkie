function [mat,x,w,u,v] = matrin(n,xs)
%LEGE.MATRIN form the matrix for interpolating from n Legendre nodes to
% the points xs
%
% optionally returns the roots, weights, and matrices for switching
% between values and coefficients (u : values -> coeffs, 
% v : coeffs -> values), which is the output of LEGE.EXPS

[x,w,u,v] = lege.exps(n);

xsvals = lege.pols(xs,n-1);

mat = (xsvals).'*u;

end