function [cf1,cf2] = besseldiff_etc_pscoefs(nterms)
% powerseries coefficients for J_0 - 1 and the other series in
% definition of Y_0
%
% these are both series with the terms z^2,z^4,z^6,..z^(2*nterms)
%
% cf1 are power series coeffs of even terms (excluding constant) 
% for J_0(z) - 1
%
% cf2 are power series coeffs for the similar series
% z^2/4 - (1+1/2)*(z^2/4)^2/(2!)^2 + (1+1/2+1/3)*(z^2/4)^3/(3!)^2 - ...
%
% which appears in the power series for H_0(z)

sum1 = 0;
fac1 = 1;
pow1 = 1;

cf1 = zeros(nterms,1);
cf2 = zeros(nterms,1);
sgn = 1;

for j = 1:nterms
    fac1 = fac1/(j*j);
    sum1 = sum1 + 1.0/j;
    pow1 = pow1*0.25;
    cf1(j) = -sgn*pow1*fac1;
    cf2(j) = sgn*sum1*pow1*fac1;
    sgn = -sgn;
end

end