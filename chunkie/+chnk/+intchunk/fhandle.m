function fint = chunkerintchunk_fhandle(f,rc,dc)

dsdtfun = @(t) sqrt(sum((lege.exev(t,dc.').^2).',1));
fintfun = @(t) f(lege.exev(t,rc.').').*dsdtfun(t);

fint = quadgk(fintfun,-1,1);

end
