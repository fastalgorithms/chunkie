function val = exev(xs,coeff)
%LEGE.EXEV evaluate legendre expansion with given coeff at points xs

sz = size(xs);
xs = xs(:);
pjm2=ones(length(xs),1);
pjm1=xs;
 
n = length(coeff);

val=coeff(1)*pjm2;
if (n == 1); return; end
val = val +coeff(2)*pjm1;
if (n == 2); return; end
 
for j = 2:(n-1)
    pj= ( (2*j-1)*xs.*pjm1-(j-1)*pjm2 ) / j;
    val = val+coeff(j+1)*pj;
    pjm2 = pjm1;
    pjm1 = pj;
end

val = reshape(val,sz);