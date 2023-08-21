function [val] = asymintdaa_v(x,r,k,efac,k2,r2)
    rdiv = r2.^2;
    prefac = 0.25*(x.^2).*efac./rdiv;
    fac = k2.^2-3*k2+3;
    val = fac.*prefac;
    %val = -exp(1i*k*sqrt(1-(1-a)*cos(x))).*sqrt(1-(1-a)*cos(x));
end

