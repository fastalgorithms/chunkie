
% test bessel0-1 and other ps coefficients 

clearvars; clc

gam = 0.57721566490153286060651209;

nterms = 12;
[cf1,cf2] = besseldiff_etc_pscoefs(nterms);

z = 1.3 + 0.1*1i;
v1 = even_pseval(cf1,z);
abs(v1-(besselj(0,z)-1))

v2 = even_pseval(cf2,z);
abs(besselh(0,z)-besselj(0,z)-1i*(2/pi)*((log(z/2)+gam)*besselj(0,z) + v2))
