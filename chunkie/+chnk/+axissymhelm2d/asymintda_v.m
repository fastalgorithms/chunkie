function [val] = asymintda_v(x,r,k,efac)
    val = x.*efac./(r.^2).*(1i*k.*r-1)/2;
    %val = exp(1i*k*sqrt(1-(1-a)*cos(x)))./sqrt(1-(1-a)*cos(x));
end

