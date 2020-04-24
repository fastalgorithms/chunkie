function val = exev(xs,coeff)
%LEGE.EXEV evaluate legendre expansion with given coeff at points xs
%
% input:
%    xs - array of points to evaluate at
%    coeff - n times nfun array of coefficients, each column corresponding
%            to an (n-1)th order Chebyshev polynomial
%
% output:
%    val - size of xs by nfun array of values for each set of coefficients
%          evaluated at the points xs (if nfun = 1, we squeeze the output)

[n,nfun] = size(coeff);
sz = size(xs);
xs = xs(:);
pjm2=ones(length(xs),1);
pjm1=xs;
 
val=pjm2*coeff(1,:);
if (n == 1); return; end
val = val + pjm1*coeff(2,:);
if (n == 2); return; end
 
for j = 2:(n-1)
    pj= ( (2*j-1)*xs.*pjm1-(j-1)*pjm2 ) / j;
    val = val+pj*coeff(j+1,:);
    pjm2 = pjm1;
    pjm1 = pj;
end

sz2 = sz;
if nfun > 1
    sz2 = [sz nfun];
end
val = squeeze(reshape(val,sz2));