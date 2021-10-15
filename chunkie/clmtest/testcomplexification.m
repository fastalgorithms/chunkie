function testcomplexification
close all
format long
format compact
% target point
xtarg = -5;
ytarg = 6;
ztarg = xtarg+1i*ytarg;
% source point
xsource = 0;
ysource = -2;
zsource = xsource + 1i*ysource;
% wave number
omega = 5.5+1i*1e-1;


% integration domain [-L,L]
L = 12;
% number of equispaced quadrature points
n = 200;
% quadrature weights
h = L/n;

% original integration domain
x = h*(-n:n);
xp = ones(size(x));

n0 = n/2;

c=12;
% contour deformation on the left
indl = 1:n0;
x0 = x(indl);
x(indl) = x0-c*1i*erfc(x0+L);
xp(indl) = 1+c*1i*exp(-(x0+L).^2)*2/sqrt(pi);

% contour deformation on the right
indr = (length(x)-n0+1):length(x);
x(indr) = -fliplr(x(indl));
xp(indr) = fliplr(xp(indl));

plot(x)
r1 = sqrt((x-xtarg).^2+ytarg^2);
% single layer potential kernel values on the contour
kerval = 1i/4*besselh(0,omega*r1);

% density
r2 = sqrt((x-xsource).^2+ysource^2);
dens = 1i/4*besselh(1,omega*r2)*omega*ysource./r2;

% trapezoidal rule 
ucomp = -2*h*sum(kerval.*dens.*xp)

% exact solution is a monopole
uexact = 1i/4*besselh(0,omega*abs(ztarg-zsource))
rerr = abs((ucomp-uexact)/uexact)

%find(imag(sqrt((x-x.').^2))<0)
