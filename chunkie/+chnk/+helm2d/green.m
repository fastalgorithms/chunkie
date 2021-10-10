function [val,grad,hess] = green(k,src,targ)
%CHNK.HELM2D.GREEN evaluate the helmholtz green's function
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
    h0 = besselh(0,1,k*r);
    val = 0.25*1i*h0;
end

[m,n] = size(xs);

if nargout > 1
%     h1 = besselh(1,1,k*r);
%     
%     grad = zeros(m,n,2);
%     
%     grad(:,:,1) = -1i*k*0.25*h1.*rx./r;
%     grad(:,:,2) = -1i*k*0.25*h1.*ry./r;
    
    h1 = besselh(1,1,k*r)./r;
    
    grad = zeros(m,n,2);
    
    ck = 0.25*1i*k;
    grad(:,:,1) = -ck*h1.*rx;
    grad(:,:,2) = -ck*h1.*ry;    
    
end

if nargout > 2

    hess = zeros(m,n,3);

%    r3 = r.^3;
    
%    h2 = 2*h1./(k*r)-h0;
%    tmp1 = (rx-ry).*(rx+ry).*h1./r3;
    h2 = 2*h1/k-h0;
    tmp1 = (rx2-ry2).*h1./r2;
    tmp2 = k*h0./r2;
    
    %hess(:,:,1) = 0.25*1i*k*((rx-ry).*(rx+ry).*h1./r3 - k*rx2.*h0./r2);
    hess(:,:,1) = ck*(tmp1 - rx2.*tmp2);
    hess(:,:,2) = ck*k*rx.*ry.*h2./r2;
    %hess(:,:,3) = 0.25*1i*k*((ry-rx).*(rx+ry).*h1./r3 - k*ry2.*h0./r2);
    hess(:,:,3) = ck*(-tmp1 - ry2.*tmp2);
end
