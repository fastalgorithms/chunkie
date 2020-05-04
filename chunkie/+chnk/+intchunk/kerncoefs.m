function fint = kerncoefs(kern,opdims,rdim,fc, ...
    rci,dci,d2ci,targ,quadgkparams)

fintfun = @(t) fkerneval(t,kern,fc,rci,dci,d2ci,targ,opdims,rdim);

fint = quadgk(fintfun,-1,1,quadgkparams{:});

end

function fvals = fkerneval(t,kern,fc,rci,dci,d2ci,targ,opdims,...
    rdim)

densvals = zeros(opdims(2),length(t));
for i = 1:opdims(2)
    densvals(i,:) = lege.exev(t,fc(:,i)); 
end

[dim,k] = size(rci);

ri = lege.exev(t,rci.').';
di = lege.exev(t,dci.').';
d2i = lege.exev(t,d2ci.').';

targinfo = []; targinfo.r = targ;
srcinfo = []; srcinfo.r = ri; srcinfo.d = di; srcinfo.d2 = d2i;
dsdt = sqrt(sum(di.^2,1));
kernmat = kern(srcinfo,targinfo);
kernmat = kernmat(rdim:opdims(1):end,:);

kernmat = reshape(kernmat,opdims(2),length(t));
fvals = kernmat.*densvals;
fvals = sum(fvals,1);
fvals = reshape(fvals,size(t)).*dsdt;

end

