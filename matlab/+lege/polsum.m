function [pol,der,tot] = polsum(xs,n)
%POLSUM evaluate nth order Legendre polynomial via recursion along with the
% sum of the squared values of P_k(xs) divided by k+1+0.5 (divided by 2 
% for P_0 and by 1.5 for P_1) --- useful for computing the roots
%
% Copyright (C) 2009: Vladimir Rokhlin
% 
% This software is being released under a modified FreeBSD license
%

% 
%        if n=0 or n=1 - exit
% 
pol=ones(size(xs));
der=zeros(size(xs));
tot = pol.^2/2.0;

if (n<=0)
    return
end

pol=xs;
der=ones(size(xs));
tot=tot+pol.^2.*(1+0.5);

if (n==1)
    return
end

pk = ones(size(xs));
pol = xs;

%       n is greater than 1. conduct recursion
% 

for k=1:(n-1)
    pkm1=pk;
    pk=pol;
    pol= ( (2*k+1)*xs.*pk-k*pkm1 )/(k+1);
    tot=tot+pol.^2.*(k+1+0.5);
end
 

%        calculate the derivative
% 

der=n*(xs.*pol-pk)./(xs.^2-1);

end
