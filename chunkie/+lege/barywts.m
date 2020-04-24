function w = barywts(k)
%LEGE.BARYWTS
%
% returns the weights for barycentric Lagrange interpolation
% for the kth order Legendre nodes
%
% input:
% 
% k - the number of Legendre nodes
%
% output:
%
% w - vector of length k such that
%
%   w(i) = 1/\prod_{j\neq i} (x(j)-x(i))
%
% where x(j) is the jth Legendre node of order k

x = lege.exps(k);

xx = bsxfun(@minus,x,x.');
xx(1:k+1:end) = 1;

w = prod(xx,1);
w = w(1)./w(:);