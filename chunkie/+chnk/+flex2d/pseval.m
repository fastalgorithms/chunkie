function [v] = pseval(cf,x)
% evaluates a power series
% with coefficients given in cf. 
%
% x can be an array

xp = ones(size(x));

v = zeros(size(x));

nterms = length(cf);
for j = 1:nterms
    v = v + cf(j)*xp;
    xp = xp.*x;
end