function [v] = even_pseval(cf,x)
% evaluates an even power series with no constant term 
% with coefficients given in cf. 
%
% x can be an array

x2 = x.*x;

xp = x2;

v = zeros(size(x));

nterms = length(cf);
for j = 1:nterms
    v = v + cf(j)*xp;
    xp = xp.*x2;
end
