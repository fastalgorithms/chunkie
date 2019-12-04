function rnorms = normals(chnkr)

assert(chnkr.dim == 2,'normals only well-defined for dim=2');
k = chnkr.k;
nch = chnkr.nch;
d = chnkr.d;

dd = sqrt(d(1,:,:).^2 + d(2,:,:).^2);
rnorms = zeros(2,k,nch);
rnorms(1,:,:) = d(2,:,:)./dd;
rnorms(2,:,:) = -d(1,:,:)./dd;

end

