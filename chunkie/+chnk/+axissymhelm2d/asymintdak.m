function [val] = asymintdak(x,a,k)
    val = -k/2*cos(x).*exp(1i*k*sqrt(1-(1-a)*cos(x)))./sqrt(1-(1-a)*cos(x));
end

