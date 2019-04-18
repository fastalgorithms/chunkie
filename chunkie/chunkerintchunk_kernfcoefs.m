function fint = chunkerintchunk_kernfcoefs(kern,ndims,rdim,fc,rci,dci,...
    targ,targtau)

fintfun = @(t) fkerneval(t,kern,fc,rci,dci,targ,targtau,ndims,rdim);

fint = quadgk(fintfun,-1,1);

end

function fvals = fkerneval(t,kern,fc,rci,dci,targ,targtau,ndims,...
    rdim)

densvals = zeros(ndims(2),length(t));
for i = 1:ndims(2)
    densvals(i,:) = lege.exev(t,fc(:,i)); 
end

[dim,k] = size(rci);
ri = zeros(dim,length(t));
di = zeros(dim,length(t));
for i = 1:dim
    ri(i,:) = lege.exev(t,rci(i,:));
    di(i,:) = lege.exev(t,dci(i,:));
end
dsdt = sqrt(sum(abs(di).^2,1));
taui = bsxfun(@rdivide,di,dsdt);

kernmat = kern(ri,targ,taui,targtau);
kernmat = kernmat(rdim:ndims(1):end,:);

kernmat = reshape(kernmat,ndims(2),length(t));
fvals = kernmat.*densvals;
fvals = sum(fvals,1);
fvals = reshape(fvals,size(t)).*dsdt;

end

