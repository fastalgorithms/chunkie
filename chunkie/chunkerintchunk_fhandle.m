function fint = chunkerinteriortegralchunk_fhandle(f,xci,yci,xpci,ypci)

dsdtfun = @(t) sqrt(lege.exev(t,xpci).^2 + lege.exev(t,ypci).^2);
fintfun = @(t) feval_handle(f,lege.exev(t,xci),lege.exev(t,yci))...
    .*dsdtfun(t);

fint = quadgk(fintfun,-1,1);

end


function feval = feval_handle(f,x,y)

xy = [ (x(:)).'; (y(:)).' ];
feval = f(xy);

end