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
    if ~isfield(targinfo,'n')
        targnorm = chnk.normal2d(targinfo);
    else
        targnorm = targinfo.n;
    end
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);

% Tangential derivative of single layer
case {'stau'}
    if ~isfield(targinfo,'n')
        targnorm = chnk.normal2d(targinfo);
    else
        targnorm = targinfo.n;
    end
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (-grad(:,:,1).*ny + grad(:,:,2).*nx);

% Hilbert transform (two times the adjoint of stau)
case {'hilb'} 
    if ~isfield(srcinfo,'n')
        srcnorm = chnk.normal2d(srcinfo);
    else
        srcnorm = srcinfo.n;
    end
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((srcnorm(1,:)),nt,1);
    ny = repmat((srcnorm(2,:)),nt,1);

    submat = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

% gradient of single layer
case {'sgrad','sg'}
    [~,grad] = chnk.lap2d.green(src,targ,true);
    submat = reshape(permute(grad,[3,1,2]),2*nt,ns);

% gradient of double layer
case {'dgrad','dg'}
    if ~isfield(srcinfo,'n')
        srcnorm = chnk.normal2d(srcinfo);
    else
        srcnorm = srcinfo.n;
    end
    [~,~,hess] = chnk.lap2d.green(src,targ,true);
    submat = -(hess(:,:,1:2).*srcnorm(1,:)+hess(:,:,2:3).*srcnorm(2,:));
    submat = reshape(permute(submat,[3,1,2]),2*nt,ns);

% normal derivative of double layer
case {'dprime','dp'}
    if ~isfield(targinfo,'n')
        targnorm = chnk.normal2d(targinfo);
    else
        targnorm = targinfo.n;
    end
    if ~isfield(srcinfo,'n')
        srcnorm = chnk.normal2d(srcinfo);
    else
        srcnorm = srcinfo.n;
    end
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
    if ~isfield(srcinfo,'n')
        srcnorm = chnk.normal2d(srcinfo);
    else
        srcnorm = srcinfo.n;
    end
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    [s,grad] = chnk.lap2d.green(src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -coef(1)*(grad(:,:,1).*nx + grad(:,:,2).*ny) + coef(2)*s;

% normal derivative of combined field
case {'cprime','cp'}
    if ~isfield(targinfo,'n')
        targnorm = chnk.normal2d(targinfo);
    else
        targnorm = targinfo.n;
    end
    if ~isfield(srcinfo,'n')
        srcnorm = chnk.normal2d(srcinfo);
    else
        srcnorm = srcinfo.n;
    end
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    [~,grad,hess] = chnk.lap2d.green(src,targ);
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submatdp = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);

    % S'
    submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
    
    submat = coef(1)*submatdp + coef(2)*submatsp;


% gradient of combined field
case {'cgrad', 'cg'}
    if ~isfield(srcinfo,'n')
        srcnorm = chnk.normal2d(srcinfo);
    else
        srcnorm = srcinfo.n;
    end
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    [~,grad,hess] = chnk.lap2d.green(src,targ,true);
    submat = -(hess(:,:,1:2).*srcnorm(1,:)+hess(:,:,2:3).*srcnorm(2,:));
    submat = coef(1)*reshape(permute(submat,[3,1,2]),2*nt,ns);
    submat = submat+coef(2)*reshape(permute(grad,[3,1,2]),2*nt,ns);



% integral of the single layer
case{'sint'}
    %srcnorm = chnk.normal2d(srcinfo);
    [s] = chnk.lap2d.green(src,targ);
    rx = bsxfun(@minus,targ(1,:).',src(1,:));
    ry = bsxfun(@minus,targ(2,:).',src(2,:));
    nrt = rx.*targinfo.n(1,:).' + ry.*targinfo.n(2,:).';
    submat = nrt.*(s/2+1/(8*pi));

% transpose of the previous one.
case{'sintt'}
    %srcnorm = chnk.normal2d(srcinfo);
    [s] = chnk.lap2d.green(src,targ);
    rx = bsxfun(@minus,targ(1,:).',src(1,:));
    ry = bsxfun(@minus,targ(2,:).',src(2,:));
    nrt = rx.*srcinfo.n(1,:) + ry.*srcinfo.n(2,:);
    submat = -nrt.*(s/2+1/(8*pi));

% integral of the double layer
case{'dint'}
    %srcnorm = chnk.normal2d(srcinfo);
    [s] = chnk.lap2d.green(src,targ);
    ntsx = bsxfun(@times,targinfo.n(1,:).',srcinfo.n(1,:));
    ntsy = bsxfun(@times,targinfo.n(2,:).',srcinfo.n(2,:));
    nts = ntsx+ntsy;
    submat = -nts.*s;
% integral of the combined field (coef(1)* integral D + coef(2)*integral S)
case{'cint'}
    coef = ones(2,1);
    if(nargin == 4); coef = varargin{1}; end
    %srcnorm = chnk.normal2d(srcinfo);
    [s] = chnk.lap2d.green(src,targ);
    ntsx = bsxfun(@times,targinfo.n(1,:).',srcinfo.n(1,:));
    ntsy = bsxfun(@times,targinfo.n(2,:).',srcinfo.n(2,:));
    nts = ntsx+ntsy;
    submat = -nts.*s;

    rx = bsxfun(@minus,targ(1,:).',src(1,:));
    ry = bsxfun(@minus,targ(2,:).',src(2,:));
    nrt = rx.*targinfo.n(1,:).' + ry.*targinfo.n(2,:).';
    submat = coef(1)*submat+coef(2)*(nrt.*(s/2+1/(8*pi)));


otherwise
    error('Unknown Laplace kernel type ''%s''.', type);
end

end
