function submat = buildmat(chnkr,kern,opdims,i,j,whts)
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
    [~,whts] = lege.exps(chnkr.k);
end

% grab specific boundary data

r = chnkr.r;
d = chnkr.d;
h = chnkr.h;

[dim,k,~] = size(r);
rs = r(:,:,j); rt = r(:,:,i); ds = d(:,:,j); dt = d(:,:,i);
rs = reshape(rs,dim,k*length(j)); rt = reshape(rt,dim,k*length(i));
ds = reshape(ds,dim,k*length(j)); dt = reshape(dt,dim,k*length(i));

hs = h(j); ht = h(i);

dsnrms = sqrt(sum(abs(ds).^2,1));
taus = bsxfun(@rdivide,ds,dsnrms);

dtnrms = sqrt(sum(abs(dt).^2,1));
taut = bsxfun(@rdivide,dt,dtnrms);

ws = kron(hs(:),whts(:));

dsdt = dsnrms(:).*ws;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

submat = bsxfun(@times,kern(rs,rt,taus,taut),(dsdtndim2(:)).');

end
