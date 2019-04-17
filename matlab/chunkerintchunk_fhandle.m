function fint = chunkerintchunk_fhandle(f,xci,yci,xpci,ypci)

dsdtfun = @(t) sqrt(legeexevvec(t,xpci).^2 + legeexevvec(t,ypci).^2);
fintfun = @(t) feval_handle(f,legeexevvec(t,xci),legeexevvec(t,yci))...
    .*dsdtfun(t);

fint = quadgk(fintfun,-1,1);

end


function feval = feval_handle(f,x,y)

xy = [ (x(:)).'; (y(:)).' ];
feval = f(xy);

end