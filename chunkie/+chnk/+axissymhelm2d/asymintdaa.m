function [val] = asymintdaa(x,a,k)
    zz = sqrt(1-(1-a)*cos(x));
    prefac = 0.25*(cos(x).^2).*exp(1i*k*zz)./(zz.^5);
    fac = -k^2*zz.^2-3*1i*k*zz+3;
    val = fac.*prefac;
    %val = -exp(1i*k*sqrt(1-(1-a)*cos(x))).*sqrt(1-(1-a)*cos(x));
end

