function [cf1,cf2] = besselkikdiff_etc_pscoefs(nterms)
% powerseries coefficients for J_0(kz)-J_0(ikz) and the other series in
% definition of Y_0(kz)-Y_0(ikz)
%
% these are both series with the terms z^2,z^6,z^10,..z^(4*nterms-2)
%
% cf1 are power series coeffs of these terms for J_0(z)-J_0(iz)
%
% cf2 are power series coeffs for the similar series
% z^2/4 + (1+1/2+1/3)*(z^2/4)^3/(3!)^2 + (1+1/2+1/3+1/4+1/5)*(z^2/4)^5/(5!)^2 + ...
%
% which appears in the power series for Y_0(z)-Y_0(iz)

sum1 = 0;
fac1 = 2;
pow1 = 1;

cf1 = zeros(nterms,1);
cf2 = zeros(nterms,1);

for j = 1:nterms
    fac1 = fac1/((2*j-1)*(2*j-1));
    sum1 = sum1 + 1.0/(2*j-1);
    pow1 = pow1*0.25;
    cf1(j) = -pow1*fac1;
    cf2(j) = sum1*pow1*fac1;
    fac1 = fac1/((2*j)*(2*j));
    sum1 = sum1 + 1.0/(2*j);
    pow1 = pow1*0.25;
end

end