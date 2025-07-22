function submat = kern(srcinfo,targinfo,type,varargin)
%CHNK.LAP2D.KERN standard Laplace layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo,targinfo,type,varargin)

src = srcinfo.r(:,:);
targ = targinfo.r(:,:);

[~,ns] = size(src);
[~,nt] = size(targ);

switch lower(type)
% double layer
case {'d', 'double'}
    %srcnorm = chnk.normal2d(srcinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    submat = -(grad(:,:,1).*srcinfo.n(1,:) + grad(:,:,2).*srcinfo.n(2,:));

% normal derivative of single layer
case {'sp', 'sprime'}
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targinfo.n(1,:)).',1,ns);
    ny = repmat((targinfo.n(2,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);

% Tangential derivative of single layer
case {'stau'}
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targinfo.n(1,:)).',1,ns);
    ny = repmat((targinfo.n(2,:)).',1,ns);

    submat = (-grad(:,:,1).*ny + grad(:,:,2).*nx);

% Hilbert transform (two times the adjoint of stau)
case {'hilb'} 
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((srcinfo.n(1,:)),nt,1);
    ny = repmat((srcinfo.n(2,:)),nt,1);

    submat = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

% gradient of single layer
case {'sgrad','sg'}
    [~,grad] = chnk.lap2d.green(src,targ,true);
    submat = reshape(permute(grad,[3,1,2]),2*nt,ns);

% gradient of double layer
case {'dgrad','dg'}
    [~,~,hess] = chnk.lap2d.green(src,targ,true);
    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = reshape(permute(submat,[3,1,2]),2*nt,ns);

% normal derivative of double layer
case {'dprime','dp'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,~,hess] = chnk.lap2d.green(src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);

% single layer
case {'s', 'single'}
    submat = chnk.lap2d.green(src,targ);

% combined field
case {'c', 'combined'}
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    [s,grad] = chnk.lap2d.green(src,targ);
    nx = repmat(srcinfo.n(1,:),nt,1);
    ny = repmat(srcinfo.n(2,:),nt,1);
    submat = -coef(1)*(grad(:,:,1).*nx + grad(:,:,2).*ny) + coef(2)*s;

% gradient of combined field
case {'cgrad', 'cg'}
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    [~,grad,hess] = chnk.lap2d.green(src,targ,true);
    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = coef(1)*reshape(permute(submat,[3,1,2]),2*nt,ns);
    submat = submat+coef(2)*reshape(permute(grad,[3,1,2]),2*nt,ns);

otherwise
    error('Unknown Laplace kernel type ''%s''.', type);
end

end

