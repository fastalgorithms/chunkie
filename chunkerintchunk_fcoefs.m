function fint = chunkerinteriortegralchunk_fcoefs(fc,xpci,ypci)

dsdtfun = @(t) sqrt(lege.exev(t,xpci).^2 + lege.exev(t,ypci).^2);
fintfun = @(t) lege.exev(t,fc).*dsdtfun(t);

fint = quadgk(fintfun,-1,1);