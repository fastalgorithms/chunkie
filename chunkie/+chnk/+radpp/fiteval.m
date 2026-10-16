function val = fiteval(r,breaks,coefs)
%CHNK.RADPP.FITEVAL evaluate a chnk.radpp.fit expansion
%
% Syntax: val = chnk.radpp.fiteval(r,breaks,coefs)
%
% breaks (nbin+1,1) and coefs (nbin,n,m) are from the info output of
% CHNK.RADPP.FIT. r is any array; val is (numel(r),m). Radii outside
% [breaks(1),breaks(end)] give NaN.
%
% see also CHNK.RADPP.FIT, CHNK.RADPP.ROTATE

[~,n,m] = size(coefs);
r = r(:);

ib   = discretize(r,breaks);
ibad = isnan(ib);
ib(ibad) = 1;

x = r - breaks(ib);
val = reshape(coefs(ib,1,:),numel(r),m);
for i = 2:n
    val = reshape(coefs(ib,i,:),numel(r),m) + x.*val;
end
val(ibad,:) = NaN;

end
