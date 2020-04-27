function fint = chunkerintchunk_fcoefs(fc,dc)


dsdtfun = @(t) sqrt(sum((lege.exev(t,dc.').^2).',1));
fintfun = @(t) lege.exev(t,fc).*dsdtfun(t);

fint = quadgk(fintfun,-1,1);