function [val] = asymintda(x,a,k)
    zz = sqrt(1-(1-a)*cos(x));
    val = cos(x).*exp(1i*k*zz)./(zz.^3).*(1i*k*zz-1)/2;
    %val = exp(1i*k*sqrt(1-(1-a)*cos(x)))./sqrt(1-(1-a)*cos(x));
end

