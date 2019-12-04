function fint = chunkerintchunk_kernfcoefs(kern,opdims,rdim,fc,rci,dci,...
    targ,targtau)

fintfun = @(t) fkerneval(t,kern,fc,rci,dci,targ,targtau,opdims,rdim);

fint = quadgk(fintfun,-1,1);

end

function fvals = fkerneval(t,kern,fc,rci,dci,targ,targtau,opdims,...
    rdim)

densvals = zeros(opdims(2),length(t));
for i = 1:opdims(2)
    densvals(i,:) = lege.exev(t,fc(:,i)); 
end

[dim,k] = size(rci);

ri = lege.exev(t,rci.').';
di = lege.exev(t,dci.').';

dsdt = sqrt(sum(abs(di).^2,1));
taui = bsxfun(@rdivide,di,dsdt);

kernmat = kern(ri,targ,taui,targtau);
kernmat = kernmat(rdim:opdims(1):end,:);

kernmat = reshape(kernmat,opdims(2),length(t));
fvals = kernmat.*densvals;
fvals = sum(fvals,1);
fvals = reshape(fvals,size(t)).*dsdt;

end

