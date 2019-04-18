
function [val,grad,hess] = glapfun(src,targ)

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

r4 = r2.*r2;

if nargout > 0

    val = -log(r2)/(4.0*pi);
    
end

[m,n] = size(xs);

if nargout > 1

    grad = zeros(m,n,2);

    grad(:,:,1) = -rx./(2.0*pi*r2);
    grad(:,:,2) = -ry./(2.0*pi*r2);
    
end

if nargout > 2

    hess = zeros(m,n,3);

    hess(:,:,1) = rx2./(pi*r4)-1.0./(2.0*pi*r2);
    hess(:,:,2) = rx.*ry./(pi*r4);
    hess(:,:,3) = ry2./(pi*r4)-1.0./(2.0*pi*r2);
end