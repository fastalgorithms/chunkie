function [pols,ders] = pols(xs,n)
%LEGE.POLS evaluate up to nth degree Legendre polynomial via recursion
%
% Syntax:
%
%   [pols,ders] = lege.pols(xs,n);
%
% Input:
%   xs - points where you want to evaluate the Legendre polynomials
%   n - highest degree polynomial of interest 
%
% Output:
%   pols - (n+1) x size(xs) array of polynomial values. pols(i,j) has the
%                 value P_{i-1}(xs(j))
%   ders - (n+1) x size(xs) array of polynomial derivative values. 
%                 ders(i,j) has the value P'_{i-1}(xs(j))
%
% see also LEGE.POL

%
% Copyright (C) 2009: Vladimir Rokhlin
% 
% This software is being released under a modified FreeBSD license
%

assert(n>=0,'n must be non-negative');
szx = size(xs);
xs = xs(:);
pols = zeros(length(xs),n+1);

pols(:,1)=ones(length(xs),1);
ders(:,1)=zeros(length(xs),1);

if (n<=0)
    
    pols = reshape(pols.',[n+1, szx]);
    ders = reshape(ders.',[n+1, szx]);
    return
end

pols(:,2)=xs(:);
ders(:,2)=ones(length(xs),1);

if (n==1)

    pols = reshape(pols.',[n+1, szx]);
    ders = reshape(ders.',[n+1, szx]);    
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
