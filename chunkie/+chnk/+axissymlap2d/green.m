function [val, grad] = green(src, targ, origin)
%CHNK.AXISSYMHELM2D.GREEN evaluate the Laplace green's function
% for the given sources and targets. 
%
% Note: that the first coordinate is r, and the second z.
% The code relies on precomputed tables and hence loops are required for 
% computing various pairwise interactions.
% Finally, the code is not efficient in the sense that val, grad, hess 
% are always internally computed independent of nargout
%
% Since the kernels are not translationally invariant in r, the size
% of the gradient is 3, for d/dr, d/dr', and d/dz
%

[~, ns] = size(src);
[~, nt] = size(targ);

gtmp = zeros(nt, ns, 3);

rt = repmat(targ(1,:).',1,ns);
rs = repmat(src(1,:),nt,1);
dz = repmat(src(2,:),nt,1)-repmat(targ(2,:).',1,ns);
r  = (rt + origin(1));
rp = (rs + origin(1));
dr = (rs-rt);
z  = zeros(size(rt));
zp = zeros(size(rt));
[g,gdz,gdr,gdrp] = chnk.axissymlap2d.gfunc(r,rp,dr,z,zp,dz);

val = g;
gtmp(:,:,1) = gdr;
gtmp(:,:,2) = gdrp;
gtmp(:,:,3) = gdz;
grad = gtmp;


