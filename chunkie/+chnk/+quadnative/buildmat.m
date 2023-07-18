function submat = buildmat(chnkr,kern,opdims,i,j,wts)
%CHNK.QUADSMOOTH.BUILDMAT build matrix for far interactions with this kernel
% assuming that the smooth rule is sufficient
% 
%


if nargin < 4
    i = 1:chnkr.nch;
end
if nargin < 5
    j = 1:chnkr.nch;
end

if nargin < 6
    [~,wts] = lege.exps(chnkr.k);
end

% grab specific boundary data

r = chnkr.rstor;
d = chnkr.dstor;
d2 = chnkr.d2stor;
h = chnkr.hstor;
n = chnkr.nstor;

[dim,k,~] = size(r);
rs = r(:,:,j); rt = r(:,:,i); ds = d(:,:,j); dt = d(:,:,i); nt = n(:,:,i);
ns = n(:,:,j);
d2s = d2(:,:,j); d2t = d2(:,:,i);
rs = reshape(rs,dim,k*length(j)); rt = reshape(rt,dim,k*length(i));
ds = reshape(ds,dim,k*length(j)); dt = reshape(dt,dim,k*length(i));
d2s = reshape(d2s,dim,k*length(j)); d2t = reshape(d2t,dim,k*length(i));

srcinfo = []; srcinfo.r = rs; srcinfo.d = ds; srcinfo.d2 = d2s; srcinfo.n = ns;
targinfo = []; targinfo.r = rt; targinfo.d = dt; targinfo.d2 = d2t; targinfo.n = nt;
hs = h(j); ht = h(i);

dsnrms = sqrt(sum(abs(ds).^2,1));
%taus = bsxfun(@rdivide,ds,dsnrms);

%dtnrms = sqrt(sum(abs(dt).^2,1));
%taut = bsxfun(@rdivide,dt,dtnrms);

ws = kron(hs(:),wts(:));

dsdt = dsnrms(:).*ws;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

submat = bsxfun(@times,kern(srcinfo,targinfo),(dsdtndim2(:)).');

end
