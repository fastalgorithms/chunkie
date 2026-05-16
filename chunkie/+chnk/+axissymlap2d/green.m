function [val, grad, hess] = green(src, targ, origin, n)
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

vtmp = zeros(nt, ns);
gtmp = zeros(nt, ns, 3);
htmp = zeros(nt, ns, 6);

rt = repmat(targ(1,:).',1,ns);
rs = repmat(src(1,:),nt,1);
dz = repmat(src(2,:),nt,1)-repmat(targ(2,:).',1,ns);
r  = (rt + origin(1));
rp = (rs + origin(1));
dr = (rs-rt);

[gval,gdz,gdr,gdrp,gdzz,gdrrp,gdrz,gdrpz] = chnk.axissymlap2d.gfunc(r,rp,dr,dz,n);

vtmp = gval;
gtmp(:,:,1) = gdr;
gtmp(:,:,2) = gdrp;
gtmp(:,:,3) = gdz;

htmp(:,:,3) = gdzz;
htmp(:,:,4) = gdrrp;
htmp(:,:,5) = gdrz;
htmp(:,:,6) = gdrpz;

if nargout > 0
    val = vtmp;
end

if nargout > 1
    grad = gtmp;
end

if nargout > 2
    hess = htmp;
end

end


