function submat = kern(srcinfo,targinfo,type,varargin)
%CHNK.LAP2D.KERN standard Laplace layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo,targinfo,type,varargin)

src = srcinfo.r(:,:);
targ = targinfo.r(:,:);

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
    %srcnorm = chnk.normal2d(srcinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    submat = -(grad(:,:,1).*srcinfo.n(1,:) + grad(:,:,2).*srcinfo.n(2,:));
end

if strcmpi(type,'sprime')
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'stau')
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (-grad(:,:,1).*ny + grad(:,:,2).*nx);
end

if strcmpi(type,'hilb') % hilbert transform (two times the adjoint of stau)
    srcnorm = chnk.normal2d(srcinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((srcnorm(1,:)),nt,1);
    ny = repmat((srcnorm(2,:)),nt,1);

    submat = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
end

if strcmpi(type,'sgrad')
    [~,grad] = chnk.lap2d.green(src,targ,true);
    submat = reshape(permute(grad,[3,1,2]),2*nt,ns);
end

if strcmpi(type,'dgrad')
    [~,~,hess] = chnk.lap2d.green(src,targ,true);
    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = reshape(permute(submat,[3,1,2]),2*nt,ns);
end

if strcmpi(type,'dprime')
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,~,hess] = chnk.lap2d.green(src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
end


if strcmpi(type,'s')
    submat = chnk.lap2d.green(src,targ);
end

if strcmpi(type,'c')
    srcnorm = chnk.normal2d(srcinfo);
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    [s,grad] = chnk.lap2d.green(src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -coef(1)*(grad(:,:,1).*nx + grad(:,:,2).*ny) + coef(2)*s;
end

if strcmpi(type,'cgrad')
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    [~,grad,hess] = chnk.lap2d.green(src,targ,true);
    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = coef(1)*reshape(permute(submat,[3,1,2]),2*nt,ns);
    submat = submat+coef(2)*reshape(permute(grad,[3,1,2]),2*nt,ns);
end
