function [xs0, ws0] = getremovablequad(k, nfac)
%CHNK.QUADGGQ.GETREMOVABLEQUAD returns a quadrature
% for integrating piecewise smooth quadrature, where the rule
% on each panel is upsampled by a factor of nfac, i.e.
% an order nfac*k quadrature rule is used on both panels

if nargin < 2
    nfac = 1;
end


kover = ceil(k*nfac);
[xleg, wleg, ~, ~] = lege.exps(k);
if kover == k
    xover = xleg;
    wover = wleg;
else
    [xover, wover, ~, ~] = lege.exps(kover);
end

xover01 = (xover+1)/2;
wover01 = wover/2;

xs0 = cell(k,1);
ws0 = cell(k,1);
for i=1:k
    xn = xleg(i);
    xls = xover01*(xn+1)-1;
    wls = wover01*(xn+1);
    xrs = xover01*(1-xn)+xn;
    wrs = wover01*(1-xn);
    xs0{i} = [xls; xrs];
    ws0{i} = [wls; wrs];
end

end
    
