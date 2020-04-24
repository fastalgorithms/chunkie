function [pol,der] = pol(xs,n)
%POL evaluate nth order Legendre polynomial via recursion
%
% Copyright (C) 2009: Vladimir Rokhlin
% 
% This software is being released under a modified FreeBSD license
%

% 
%        if n=0 or n=1 - exit
% 

if (n<=0)
    pol=ones(size(xs));
    der=zeros(size(xs));
    return
end

if (n==1)
    pol=xs;
    der=ones(size(xs));
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
end
 

%        calculate the derivative
% 

der=n*(xs.*pol-pk)./(xs.^2-1);

end
