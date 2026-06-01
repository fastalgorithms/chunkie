function [val,grad,hess] = green(zkE,src,targ)
%CHNK.HELM1D.GREEN evaluate the 1D helmholtz green's function
% for the given sources and targets

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

rx = xt-xs;
ry = yt-ys;

rx2 = rx.*rx;
ry2 = ry.*ry;

r2 = rx2+ry2;

r = sqrt(r2);


if nargout > 0
    val = exp(1j*zkE*r);
end

[m,n] = size(xs);

if nargout > 1
    grad = zeros(m,n,2);
    
    grad(:,:,1) = (1/2).*(rx./r).*exp(1j*zkE*r);
    grad(:,:,2) = (1/2).*(ry./r).*exp(1j*zkE*r);
    grad = 2*1j*zkE*grad;
end

if nargout > 2

    hess = zeros(m,n,3);

    r3 = r.^3;
    
    hess(:,:,1) = (1./(2*r3)).*(1j*zkE*r.*rx2 - rx2 + r2).*exp(1j*zkE*r);
    hess(:,:,2) = (1./(2*r3)).*(1j*zkE*r.*rx.*ry - rx.*ry).*exp(1j*zkE*r);
    hess(:,:,3) = (1./(2*r3)).*(1j*zkE*r.*ry2 - ry2 + r2).*exp(1j*zkE*r);
    hess = 2*1j*zkE*hess;
end
