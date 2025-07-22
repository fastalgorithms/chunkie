function submat= kern(zk,srcinfo,targinfo,type,varargin)
%CHNK.HELM2D.KERN standard Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnk.heml2d.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) = i/4 H_0^{(1)}(zk |x-y|)
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
%                type == 'stau', tangential derivative of single layer
%                type == 'all', returns all four layer potentials, 
%                       [coef(1,1)*D coef(1,2)*S; coef(2,1)*D' coef(2,2)*S']
%                type == 'c2trans' returns the combined field, and the 
%                          normal derivative of the combined field
%                        [coef(1)*D + coef(2)*S; coef(1)*D' + coef(2)*S']
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
%                append '_diff' to type to compute difference stably
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
% see also CHNK.HELM2D.GREEN
  
src = srcinfo.r(:,:);
targ = targinfo.r(:,:);

[~,ns] = size(src);
[~,nt] = size(targ);

switch lower(type)
% double layer
case {'d', 'double'}
  srcnorm = srcinfo.n(:,:);
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

% double layer (difference)
case {'d_diff', 'double_diff'}
  srcnorm = srcinfo.n(:,:);
  [~,grad] = chnk.helm2d.helmdiffgreen(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

% normal derivative of single layer
case {'sp', 'sprime'}
  targnorm = targinfo.n(:,:);
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);

% normal derivative of single layer (difference)
case {'sp_diff', 'sprime_diff'}
  targnorm = targinfo.n(:,:);
  [~,grad] = chnk.helm2d.helmdiffgreen(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);


% Tangential derivative of single layer
case {'stau', 'st'}
  targtan = targinfo.d(:,:);
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  dx = repmat((targtan(1,:)).',1,ns);
  dy = repmat((targtan(2,:)).',1,ns);
  ds = sqrt(dx.*dx+dy.*dy);
  submat = (grad(:,:,1).*dx + grad(:,:,2).*dy)./ds;

% Tangential derivative of single layer (difference)
case {'stau_diff','st_diff'}
  targtan = targinfo.d(:,:);
  [~,grad] = chnk.helm2d.helmdiffgreen(zk,src,targ);
  dx = repmat((targtan(1,:)).',1,ns);
  dy = repmat((targtan(2,:)).',1,ns);
  ds = sqrt(dx.*dx+dy.*dy);
  submat = (grad(:,:,1).*dx + grad(:,:,2).*dy)./ds;

% single layer
case {'s', 'single'}
  submat = chnk.helm2d.green(zk,src,targ);

% single layer (difference)
case {'s_diff', 'single_diff'}
  submat = chnk.helm2d.helmdiffgreen(zk,src,targ);

% gradient of single layer
case {'sgrad','sg'}
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    submat = reshape(permute(grad,[3,1,2]),[],ns);

% gradient of single layer(difference)
case {'sgrad_diff','sg_diff'}
    [~,grad] = chnk.helm2d.helmdiffgreen(zk,src,targ);
    submat = reshape(permute(grad,[3,1,2]),[],ns);

% normal derivative of double layer
case {'dp', 'dprime'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,~,hess] = chnk.helm2d.green(zk,src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);

% normal derivative of double layer (difference)
case {'dp_diff', 'dprime_diff'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,~,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);

% gradient of double layer
case {'dgrad','dg'}
    [~,~,hess] = chnk.helm2d.green(zk,src,targ);
    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = reshape(permute(submat,[3,1,2]),2*nt,ns);

% gradient of double layer (difference)
case {'dgrad_diff','dg_diff'}
    [~,~,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);
    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = reshape(permute(submat,[3,1,2]),2*nt,ns);

% Combined field 
case {'c', 'combined'}
  srcnorm = srcinfo.n(:,:);
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*submats;

% Combined field (difference)
case {'c_diff', 'combined_diff'}
  srcnorm = srcinfo.n(:,:);
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  [submats,grad] = chnk.helm2d.helmdiffgreen(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*submats;

% normal derivative of combined field
case {'cp', 'cprime'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [~,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);

  submat = coef(1)*submatdp + coef(2)*submatsp;

% normal derivative of combined field (difference)
case {'cp_diff', 'cprime_diff'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [~,grad,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);

  submat = coef(1)*submatdp + coef(2)*submatsp;

% gradient of combined field
case {'cgrad','cg'}
    coef = ones(2,1);
    if(nargin == 5); coef = varargin{1}; end
    [~,grad,hess] = chnk.helm2d.green(zk,src,targ);

    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = coef(1)*reshape(permute(submat,[3,1,2]),2*nt,ns);
    submat = submat+coef(2)*reshape(permute(grad,[3,1,2]),2*nt,ns);

% gradient of combined field (difference)
case {'cgrad_diff','cg_diff'}
    coef = ones(2,1);
    if(nargin == 5); coef = varargin{1}; end
    [~,grad,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);

    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = coef(1)*reshape(permute(submat,[3,1,2]),2*nt,ns);
    submat = submat+coef(2)*reshape(permute(grad,[3,1,2]),2*nt,ns);

% Dirichlet and Neumann data corresponding to combined field
case {'c2trans', 'c2t', 'c2tr'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(2*nt,ns);

  submat(1:2:2*nt,:) = coef(1)*submatd + coef(2)*submats;
  submat(2:2:2*nt,:) = coef(1)*submatdp + coef(2)*submatsp;


% Dirichlet and Neumann data corresponding to combined field (difference)
case {'c2trans_diff', 'c2t_diff', 'c2tr_diff'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(2*nt,ns);

  submat(1:2:2*nt,:) = coef(1)*submatd + coef(2)*submats;
  submat(2:2:2*nt,:) = coef(1)*submatdp + coef(2)*submatsp;

% all kernels, [c11 D, c12 S; c21 D', c22 S'] 
case {'all','trans_sys','ts'}
 
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  cc = varargin{1};
  
  submat = zeros(2*nt,2*ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  submat(1:2:2*nt,1:2:2*ns) = submatd*cc(1,1);
  submat(1:2:2*nt,2:2:2*ns) = submats*cc(1,2);
  submat(2:2:2*nt,1:2:2*ns) = submatdp*cc(2,1);
  submat(2:2:2*nt,2:2:2*ns) = submatsp*cc(2,2);

% all kernels, [c11 D, c12 S; c21 D', c22 S'] (difference)
case {'all_diff','trans_sys_diff','ts_diff'}
 
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  cc = varargin{1};
  
  submat = zeros(2*nt,2*ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  
  submat(1:2:2*nt,1:2:2*ns) = submatd*cc(1,1);
  submat(1:2:2*nt,2:2:2*ns) = submats*cc(1,2);
  submat(2:2:2*nt,1:2:2*ns) = submatdp*cc(2,1);
  submat(2:2:2*nt,2:2:2*ns) = submatsp*cc(2,2);

% Dirichlet data/potential correpsonding to transmission rep
case {'trans_rep','tr'} 

  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  srcnorm = srcinfo.n(:,:);
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

  submat = zeros(1*nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatd;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submats;

% Dirichlet data/potential correpsonding to transmission rep (difference)
case {'trans_rep_diff','tr_diff'} 

  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  srcnorm = srcinfo.n(:,:);
  [submats,grad] = chnk.helm2d.helmdiffgreen(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

  submat = zeros(1*nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatd;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submats;

% Neumann data corresponding to transmission rep
case {'trans_rep_prime','trp', 'trans_rep_p'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  
  submat = zeros(nt,ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submat = zeros(nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatdp;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submatsp;

% Neumann data corresponding to transmission rep (difference)
case {'trans_rep_prime_diff','trp_diff', 'trans_rep_p_diff'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  
  submat = zeros(nt,ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submat = zeros(nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatdp;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submatsp;

% Gradient correpsonding to transmission rep
case {'trans_rep_grad','trg','trans_rep_g'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  
  srcnorm = srcinfo.n(:,:);
  
  submat = zeros(nt,ns,6);
  % S
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,1).*nxsrc + hess(:,:,2).*nysrc);
  submat(1:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,1);
    
  submat(2:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,2).*nxsrc + hess(:,:,3).*nysrc);
  submat(2:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,2);


% Gradient correpsonding to transmission rep (difference)
case {'trans_rep_grad_diff','trg_diff','trans_rep_g_diff'}
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  
  srcnorm = srcinfo.n(:,:);
  
  submat = zeros(nt,ns,6);
  % S
  [submats,grad,hess] = chnk.helm2d.helmdiffgreen(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,1).*nxsrc + hess(:,:,2).*nysrc);
  submat(1:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,1);
    
  submat(2:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,2).*nxsrc + hess(:,:,3).*nysrc);
  submat(2:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,2);
otherwise
    error('Unknown Helmholtz kernel type ''%s''.', type);
end