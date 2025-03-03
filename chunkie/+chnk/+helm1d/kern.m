function submat= kern(zkE,srcinfo,targinfo,type,varargin)
%CHNK.HELM1D.KERN standard Helmholtz layer potential kernels in 1D
% 
% Syntax: submat = chnk.heml1d.kern(zkE,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) = e^{iE|x-y|}
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
%
% Input:
%   zkE - complex number, Helmholtz wave number
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
% see also CHNK.HELM1D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

% double layer
if strcmpi(type,'d')
  srcnorm = srcinfo.n;
  [~,grad] = chnk.helm1d.green(zkE,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

% normal derivative of single layer
if strcmpi(type,'sprime')
  targnorm = targinfo.n;
  [~,grad] = chnk.helm1d.green(zkE,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end


% Tangential derivative of single layer
if strcmpi(type,'stau')
  targtan = targinfo.d;
  [~,grad] = chnk.helm1d.green(zkE,src,targ);
  dx = repmat((targtan(1,:)).',1,ns);
  dy = repmat((targtan(2,:)).',1,ns);
  ds = sqrt(dx.*dx+dy.*dy);
  submat = (grad(:,:,1).*dx + grad(:,:,2).*dy)./ds;
end

% single layer
if strcmpi(type,'s')
  submat = chnk.helm1d.green(zkE,src,targ);
end

% normal derivative of double layer
if strcmpi(type,'dprime')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  [~,~,hess] = chnk.helm1d.green(zkE,src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
end

% Combined field 
if strcmpi(type,'c')
  srcnorm = srcinfo.n;
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  [submats,grad] = chnk.helm1d.green(zkE,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*submats;
end

% normal derivative of combined field
if strcmpi(type,'cprime')
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  


  % Get gradient and hessian info
  [~,grad,hess] = chnk.helm1d.green(zkE,src,targ);

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
end


% Dirichlet and neumann data corresponding to combined field
if strcmpi(type,'c2trans') 
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm1d.green(zkE,src,targ);

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

end


% all kernels, [c11 D, c12 S; c21 D', c22 S'] 
if strcmpi(type,'all')
 
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  cc = varargin{1};
  
  submat = zeros(2*nt,2*ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm1d.green(zkE,src,targ);

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
end

% Dirichlet data/potential correpsonding to transmission rep
if strcmpi(type,'trans_rep') 

  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  srcnorm = srcinfo.n;
  [submats,grad] = chnk.helm1d.green(zkE,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

  submat = zeros(1*nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatd;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submats;
end

% Neumann data corresponding to transmission rep
if strcmpi(type,'trans_rep_prime')
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm1d.green(zkE,src,targ);

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
  submat = zeros(nt,2*ns)
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatdp;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submatsp;
end


% Gradient correpsonding to transmission rep
if strcmpi(type,'trans_rep_grad')
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,ns,6);
  % S
  [submats,grad,hess] = chnk.helm1d.green(zkE,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,1).*nxsrc + hess(:,:,2).*nysrc);
  submat(1:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,1);
    
  submat(2:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,2).*nxsrc + hess(:,:,3).*nysrc);
  submat(2:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,2);
end


