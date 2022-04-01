function submat = kern(srcinfo,targinfo,type,varargin)
%CHNK.LAP2D.KERN standard Laplace layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo,targinfo,type,varargin)

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
    %srcnorm = chnk.normal2d(srcinfo);
    srcnorm = srcinfo.n;
    [~,grad] = chnk.lap2d.green(src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sprime')
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'s')
    submat = chnk.lap2d.green(src,targ);
end

if strcmpi(type,'c')
    srcnorm = chnk.normal2d(srcinfo);
    eta = varargin{1};
    [s,grad] = chnk.lap2d.green(src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny) + eta*s;
end
