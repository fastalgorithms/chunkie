function dercoeffs = derpol(coeffs)
%LEGE.DERPOL compute the coefficients of the derivative
% of the polynomial with the coefficients given on 
% input

sz = size(coeffs);
sz(1)= max(sz(1)-1,0);
n = sz(1);
dercoeffs = zeros(sz);

if (n<=0); return; end

pk = coeffs(n+1,:); 
pkm1 = coeffs(end-1,:);
pkm2=0;

for k = n+1:-1:2
    j = k-1;
    dercoeffs(k-1,:)= pk*(2*j-1);
    if(k > 2); pkm2=coeffs(k-2,:)+pk; end
    pk=pkm1;
    pkm1=pkm2;
end

end