function [rk0,re0] = gaus_agm(x)

eps  = 1E-15;
a    = sqrt(2./(x+1));
delt = 1./sqrt(1-a.*a);
aa0  = delt + sqrt(delt.*delt-1);
bb0   = 1./(delt+sqrt(delt.*delt-1));
a0  = ones(size(delt));
b0  = 1./delt;


fact = ((a0+b0)/2).^2;

for i=1:1000
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);

    aa1 = (aa0+bb0)/2;
    bb1 = sqrt(aa0.*bb0);
    a0 = a1;
    b0 = b1;
    aa0 = aa1;
    bb0 = bb1;

    c0 = (a1-b1)/2;
    fact = fact-(c0.*c0)*2^(i);
    drel = abs(a0-b0)./abs(a0);
    drel2= abs(aa0-bb0)./abs(aa0);
    if (max(drel+drel2) <2*eps)
        break
    end

end

rk0 = pi./(2*aa0.*sqrt(1-a.*a));
re0 = pi*fact./(2*a0);

end