function zell = bernstein_ellipse(ntheta,rho)

th = linspace(0,2*pi,ntheta+1);
zell = rho*exp(1i*th(1:end-1));
zell = (zell+1./zell)/2;
zell = zell(:);