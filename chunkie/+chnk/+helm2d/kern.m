
function submat = kern(zk,srcinfo,targinfo,type,varargin)
%CHNK.HELM2D.KERN standard Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnkr.heml2d.kern(zk,srcinfo,targingo,type,varargin)
% 

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
    srcnorm = chnk.normal2d(srcinfo);
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sprime')
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'s')
    submat = chnk.helm2d.green(zk,src,targ);
end

if strcmpi(type,'c')
    srcnorm = chnk.normal2d(srcinfo);
    eta = varargin{1};
    [s,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny) + 1i*eta*s;
end
