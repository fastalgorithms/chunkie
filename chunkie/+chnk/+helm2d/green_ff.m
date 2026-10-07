function [val,grad,hess] = green_ff(k,src,targ)
%CHNK.HELM2D.GREEN_FF evaluate the helmholtz far-field green's function
% for the given sources and targets

prefac = 1i/(2*sqrt(2*pi))*exp(-1i*pi/4);

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);
rt = sqrt(xt.^2+yt.^2);

xt = xt./rt;
yt = yt./rt;

dp = xt.*xs+yt.*ys;

g0 = prefac*exp(-1i*dp*k);


if nargout > 0
    val = g0;  
end
if nargout > 1
    grad(:,:,1) = -1i*k*g0.*xt;
    grad(:,:,2) = -1i*k*g0.*yt;
    grad = -grad;
end
if nargout > 2
    hess(:,:,1) = (-1i*k)^2*g0.*xt.^2;
    hess(:,:,2) = (-1i*k)^2*g0.*xt.*yt;
    hess(:,:,3) = (-1i*k)^2*g0.*yt.^2;
end