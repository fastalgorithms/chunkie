
function [val,grad,hess] = green(zk,src,targ)

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
    h0 = besselh(0,1,zk*r);
    val = 0.25*1i*h0;
end

[m,n] = size(xs);

if nargout > 1
    h1 = besselh(1,1,zk*r);
    grad = zeros(m,n,2);
    
    grad(:,:,1) = -1i*zk*0.25*h1.*rx./r;
    grad(:,:,2) = -1i*zk*0.25*h1.*ry./r;
    
end

if nargout > 2

    hess = zeros(m,n,3);

    r3 = r.^3;
    
    h2 = 2*h1./(zk*r)-h0;
    
    hess(:,:,1) = 0.25*1i*zk*((rx-ry).*(rx+ry).*h1./r3 - zk*rx2.*h0./r2);
    hess(:,:,2) = 0.25*1i*zk*zk*rx.*ry.*h2./r2;
    hess(:,:,3) = 0.25*1i*zk*((ry-rx).*(rx+ry).*h1./r3 - zk*ry2.*h0./r2);
end
