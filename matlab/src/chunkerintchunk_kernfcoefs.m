function fint = chunkerintchunk_kernfcoefs(kern,ndims,dim,fc,xci,yci, ...
    xpci,ypci,targ,targn)

fintfun = @(t) fkerneval(t,kern,fc,xci,yci,xpci,ypci,targ,targn,ndims,dim);

fint = quadgk(fintfun,-1,1);

end

function fvals = fkerneval(t,kern,fc,xci,yci,xpci,ypci,targ,targn,ndims,...
    dim)

densvals = zeros(ndims(2),length(t));
for i = 1:ndims(2)
    densvals(i,:) = legeexevvec(t,fc(:,i)); 
end
x = legeexevvec(t,xci); x = (x(:)).';
xp = legeexevvec(t,xpci); xp = (xp(:)).';
y = legeexevvec(t,yci); y = (y(:)).';
yp = legeexevvec(t,ypci); yp = (yp(:)).';
dsdt = sqrt(xp.^2+yp.^2);
s = [ x; y ];
sn = [ yp./dsdt; -xp./dsdt ];

kernmat = kern(s,targ,sn,targn);
kernmat = kernmat(dim:ndims(1):end,:);

kernmat = reshape(kernmat,ndims(2),length(t));
fvals = kernmat.*densvals;
fvals = sum(fvals,1);
fvals = reshape(fvals,size(t)).*dsdt;

end

