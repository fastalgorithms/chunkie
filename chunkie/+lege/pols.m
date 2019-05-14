function [pols,ders] = pols(xs,n)
%LEGE.POLS evaluate up to nth order Legendre polynomial via recursion
%
% output is size(xs) x (n+1) array
%
% Copyright (C) 2009: Vladimir Rokhlin
% 
% This software is being released under a modified FreeBSD license
%

% 
%        if n=0 or n=1 - exit
% 

assert(n>=0,'n must be non-negative');
szx = size(xs);
xs = xs(:);
pols = zeros(length(xs),n+1);

pols(:,1)=ones(length(xs),1);
ders(:,1)=zeros(length(xs),1);

if (n<=0)
    return
end

pols(:,2)=xs(:);
ders(:,2)=ones(length(xs),1);

if (n==1)
    return
end

pk = ones(length(xs),1);

%       n is greater than 1. conduct recursion
% 

xs2m1 = xs.^2-1.0;
pkp1 = pols(:,2);

for k=1:(n-1)
    pkm1=pk;
    pk = pkp1;
    pkp1 = ( (2*k+1)*xs.*pk-k*pkm1 )/(k+1);
    pols(:,k+2) = pkp1;
    ders(:,k+2) = (k+1)*(xs.*pkp1-pk)./xs2m1;
end
 
pols = reshape(pols.',[n+1, szx]);
ders = reshape(ders.',[n+1, szx]);

end
