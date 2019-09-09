function rnorms = normals(trap)

assert(trap.dim == 2,'normals only well-defined for dim=2');

d = trap.d;

dd = sqrt(d(1,:,:).^2 + d(2,:,:).^2);
rnorms = zeros(2,trap.npt);
rnorms(1,:) = d(2,:)./dd;
rnorms(2,:) = -d(1,:)./dd;

end

