function [val] = asymintdk(x,a,k)
    val = 1i*exp(1i*k*sqrt(1-(1-a)*cos(x)));
end

