function submat=kern(zk,srcinfo,targinfo,type,quas_param,varargin)
%CHNK.HELM2DQUAS.KERN quasi-periodic Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnk.helm2dquas.kern(zk,srcinfo,targinfo,type,quas_param,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). 
%  
% Kernels based on the quasi periodic green's function satisfying 
%       G(x+d e_1,y) = G(x,y) exp(i kapppa d),
%  which is given by
%       G(x,y) = sum_{n=-inf}^inf i/4 H_0^{(1)}(zk |x- n d e_1-y|) 
%                   exp(i kapppa d n)
%
% NOTE: The quasiperiodic green's has extra singularities at periodic 
% copies of the source, which might cause quadrature errors.
%
% NOTE: This will be singular if kappa = +- pi/d
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'dprime', normal derivative of double layer D'
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
%                type == 'all', returns all four layer potentials, 
%                       [coef(1,1)*D coef(1,2)*S; coef(2,1)*D' coef(2,2)*S']
%                type == 'c2trans' returns the combined field, and the 
%                          normal derivative of the combined field
%                        [coef(1,1)*D + coef(1,2)*S; coef(2,1)*D' + coef(2,2)*S']
%                type == 'trans_rep' returns the potential corresponding
%                           to the transmission representation
%                        [coef(1)*D coef(2)*S]
%                type == 'trans_rep_prime' returns the normal derivative
%                          corresponding to the transmission representation
%                        [coef(1)*D' coef(2)*S']
%                type == 'trans_rep_grad' returns the gradient corresponding
%                         to the transmission representation
%                        [coef(1)*d_x D coef(2)*d_x S;
%                         coef(1)*d_y D coef(2)*d_y S]
% 
%   quas_param - struct with quasiperiodic parameters,
%       quas_param.kappa - phase differences
%       quas_param.d - period
%       quas_param.l - radius excluded from local quasiperiodic farfield
%       quas_param.sn - precomputed lattice sum integrals, see 
%               chnk.helm2dquas.latticecoefs
%
%   varargin{1} - coef: length 2 array in the combined layer 
%                 formula, 2x2 matrix for all kernels
%                 otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2DQUAS.GREEN

src = srcinfo.r(:,:);
targ = targinfo.r(:,:);


kappa = quas_param.kappa;
d = quas_param.d;
l = quas_param.l;
sn = quas_param.sn;


[~,ns] = size(src);
[~,nt] = size(targ);
nkappa = length(kappa);

% double layer
if strcmpi(type,'d')
  srcnorm = srcinfo.n(:,:);
  [~,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);
  nx = repmat(srcnorm(1,:),nkappa*nt,1);
  ny = repmat(srcnorm(2,:),nkappa*nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

% normal derivative of single layer
if strcmpi(type,'sprime')
  targnorm = targinfo.n(:,:);
  [~,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);
  nx = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
  nx = reshape(nx,[],ns);
  ny = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
  ny = reshape(ny,[],ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

% single later
if strcmpi(type,'s')
  submat = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);
  submat = reshape(submat,[],ns);
end

% normal derivative of double layer
if strcmpi(type,'dprime')
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,~,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);
  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  % nxtarg = repmat((targnorm(1,:)).',1,ns);
  % nytarg = repmat((targnorm(2,:)).',1,ns);
  nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
  nxtarg = reshape(nxtarg,[],ns);
  nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
  nytarg = reshape(nytarg,[],ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
end


% gradient of single layer
if strcmpi(type,'sgrad')
    [~,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);
    submat = reshape(permute(grad,[1,3,2]),[],ns);
end

% gradient of double layer
if strcmpi(type,'dgrad')
    [~,~,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);
    submat = -(hess(:,:,1:2).*reshape(srcinfo.n(1,:),1,[],1)+hess(:,:,2:3).*reshape(srcinfo.n(2,:),1,[],1));
    submat = reshape(permute(submat,[1,3,2]),[],ns);
end

% Combined field 
if strcmpi(type,'c')
  srcnorm = srcinfo.n(:,:);
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end
  [submats,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);
  nx = repmat(srcnorm(1,:),nkappa*nt,1);
  ny = repmat(srcnorm(2,:),nkappa*nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*reshape(submats,[],ns);
end

% normal derivative of combined field
if strcmpi(type,'cprime')
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [~,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);

  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
  nxtarg = reshape(nxtarg,[],ns);
  nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
  nytarg = reshape(nytarg,[],ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);

  submat = coef(1)*submatdp + coef(2)*submatsp;
end

% Dirichlet and neumann data corresponding to combined field
if strcmpi(type,'c2trans') 
  coef = ones(1,2);
  if(nargin >= 6); coef = varargin{1}; end
  if numel(coef) == 2, coef = repmat(coef(:).',2,1); end
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);

  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
  nxtarg = reshape(nxtarg,[],ns);
  nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
  nytarg = reshape(nytarg,[],ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  % submat = zeros(nkappa*2*nt,ns);
  % submat(1:2:end,:) = coef(1)*submatd + coef(2)*reshape(submats,[],ns);
  % submat(2:2:end,:) = coef(1)*submatdp + coef(2)*submatsp;

  submat = zeros(nkappa,2,nt,ns);
  submat(:,1,:,:) = reshape(coef(1,1)*submatd + coef(1,2)*reshape(submats,[],ns),nkappa,1,nt,[]);
  submat(:,2,:,:) = reshape(coef(2,1)*submatdp + coef(2,2)*submatsp,nkappa,1,nt,[]);
  submat = reshape(submat, [],ns);
end

% gradient of combined field
if strcmpi(type,'cgrad')
    coef = ones(2,1);
    if(nargin >= 6); coef = varargin{1}; end
    [~,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);

    submats = reshape(permute(grad,[1,3,2]),[],ns);

    submatd = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submatd = reshape(permute(submatd,[1,3,2]),[],ns);
    
    submat = coef(1)*submatd + coef(2)*submats;
end
  

% all kernels, [c11 D, c12 S; c21 D', c22 S'] 
if strcmpi(type,'all')
 
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  cc = varargin{1};

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);

  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
  nxtarg = reshape(nxtarg,[],ns);
  nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
  nytarg = reshape(nytarg,[],ns);

  submats = reshape(submats,[],ns);

  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  submat = zeros(nkappa,2,nt,2*ns);
  submat(:,1,:,1:2:2*ns) = reshape(submatd,nkappa,1,nt,[])*cc(1,1);
  submat(:,1,:,2:2:2*ns) = reshape(submats,nkappa,1,nt,[])*cc(1,2);
  submat(:,2,:,1:2:2*ns) = reshape(submatdp,nkappa,1,nt,[])*cc(2,1);
  submat(:,2,:,2:2:2*ns) = reshape(submatsp,nkappa,1,nt,[])*cc(2,2);
  submat = reshape(submat, [],2*ns);
end
% Dirichlet data/potential correpsonding to transmission rep
if strcmpi(type,'trans_rep') 

  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end;
  srcnorm = srcinfo.n(:,:);
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);

  submats = reshape(submats,[],ns);
  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  submatd = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(nkappa*nt,2*ns,'like',1i);
  submat(1:1:end,1:2:2*ns) = coef(1)*submatd;
  submat(1:1:end,2:2:2*ns) = coef(2)*submats;
end

% Neumann data corresponding to transmission rep
if strcmpi(type,'trans_rep_prime')
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end;
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  
  submat = zeros(nt,ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);

  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
  nxtarg = reshape(nxtarg,[],ns);
  nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
  nytarg = reshape(nytarg,[],ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submat = zeros(nkappa*nt,2*ns);
  submat(1:1:end,1:2:2*ns) = coef(1)*submatdp;
  submat(1:1:end,2:2:2*ns) = coef(2)*submatsp;
end

% Gradient correpsonding to transmission rep
if strcmpi(type,'trans_rep_grad')
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end;
  
  srcnorm = srcinfo.n(:,:);
  
  % submat = zeros(nt,ns,6);
  % S
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l);

  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(nkappa,2*nt,2*ns);
  
  submat(:,1,1:2:2*ns) = -coef(1)*(hess(:,:,1).*nxsrc + hess(:,:,2).*nysrc);
  submat(:,1,2:2:2*ns) = coef(2)*grad(:,:,1);
    
  submat(:,2,1:2:2*ns) = -coef(1)*(hess(:,:,2).*nxsrc + hess(:,:,3).*nysrc);
  submat(:,2,2:2:2*ns) = coef(2)*grad(:,:,2);
  submat = reshape(submat,nkappa*2*nt,2*ns);
end






end