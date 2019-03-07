function fint = chunkerintchunk_fcoefs(fc,xpci,ypci)

dsdtfun = @(t) sqrt(legeexevvec(t,xpci).^2 + legeexevvec(t,ypci).^2);
fintfun = @(t) legeexevvec(t,fc).*dsdtfun(t);

fint = quadgk(fintfun,-1,1);