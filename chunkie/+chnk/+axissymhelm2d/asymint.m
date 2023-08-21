function [val] = asymint(x,a,k)
    val = exp(1i*k*sqrt(1-(1-a).*cos(x)))./sqrt(1-(1-a).*cos(x));
end

