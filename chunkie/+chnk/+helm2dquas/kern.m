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
%   varargin{1} - ising: if set to 0, only include the periodic copies. 
%                 If set to 1, include the free-space part
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


ising = 1;
if length(varargin) >1
    ising = varargin{2};
end 

kappa = quas_param.kappa;
d = quas_param.d;
l = quas_param.l;
sn = quas_param.sn;


[~,ns] = size(src);
[~,nt] = size(targ);
nkappa = length(kappa);

switch lower(type)
% double layer
case {'d', 'double'}
  srcnorm = srcinfo.n(:,:);
  [~,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);
  nx = repmat(srcnorm(1,:),nkappa*nt,1);
  ny = repmat(srcnorm(2,:),nkappa*nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

% normal derivative of single layer
case {'sp', 'sprime'}
  targnorm = targinfo.n(:,:);
  [~,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);
  nx = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
  nx = reshape(nx,[],ns);
  ny = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
  ny = reshape(ny,[],ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);



% single later
case {'s', 'single'}
  submat = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);
  submat = reshape(submat,[],ns);

% normal derivative of double layer
case {'dp', 'dprime'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,~,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);
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


% gradient of single layer
case {'sgrad','sg'}
    [~,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);
    submat = reshape(permute(grad,[1,3,2]),[],ns);

% gradient of double layer
case {'dgrad','dg'}
    [~,~,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);
    submat = -(hess(:,:,1:2).*reshape(srcinfo.n(1,:),1,[],1)+hess(:,:,2:3).*reshape(srcinfo.n(2,:),1,[],1));
    submat = reshape(permute(submat,[1,3,2]),[],ns);

% Combined field 
case {'c', 'combined'}
  srcnorm = srcinfo.n(:,:);
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end
  [submats,grad] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);
  nx = repmat(srcnorm(1,:),nkappa*nt,1);
  ny = repmat(srcnorm(2,:),nkappa*nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*reshape(submats,[],ns);

% normal derivative of combined field
case {'cp', 'cprime'}
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [~,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);

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


% Dirichlet and Neumann data corresponding to combined field
case {'c2trans', 'c2t', 'c2tr'}
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end
  if numel(coef) == 2, coef = repmat(coef(:).',2,1); end

  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);

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

% gradient of combined field
case {'cgrad','cg'}
    coef = ones(2,1);
    if(nargin >= 6); coef = varargin{1}; end
    [~,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);

    submats = reshape(permute(grad,[1,3,2]),[],ns);

    submatd = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submatd = reshape(permute(submatd,[1,3,2]),[],ns);
    
    submat = coef(1)*submatd + coef(2)*submats;
  

% all kernels, [c11 D, c12 S; c21 D', c22 S'] 
case {'all','trans_sys','tsys'}
 
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  cc = varargin{1};

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);

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
% end
% % Dirichlet data/potential correpsonding to transmission rep
% if strcmpi(type,'trans_rep') 
% =======
%   submat = zeros(nkappa*2*nt,2*ns);
%   submat(1:2:end,1:2:2*ns) = submatd*cc(1,1);
%   submat(1:2:end,2:2:2*ns) = submats*cc(1,2);
%   submat(2:2:end,1:2:2*ns) = submatdp*cc(2,1);
%   submat(2:2:end,2:2:2*ns) = submatsp*cc(2,2);
% >>>>>>> master

% Dirichlet data/potential correpsonding to transmission rep
case {'trans_rep','trep'} 
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end;
  srcnorm = srcinfo.n(:,:);
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);

  submats = reshape(submats,[],ns);
  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
  submatd = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(nkappa*nt,2*ns,'like',1i);
  submat(1:1:end,1:2:2*ns) = coef(1)*submatd;
  submat(1:1:end,2:2:2*ns) = coef(2)*submats;

% Neumann data corresponding to transmission rep
case {'trans_rep_prime','trep_p', 'trans_rep_p'}
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end;
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  
  submat = zeros(nt,ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);

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

% Gradient correpsonding to transmission rep
case {'trans_rep_grad','trep_g', 'trans_rep_g'}
  coef = ones(2,1);
  if(nargin >= 6); coef = varargin{1}; end;
  
  srcnorm = srcinfo.n(:,:);

  [submats,grad,hess] = chnk.helm2dquas.green(src,targ,zk,kappa,d,sn,l,ising);

  nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
  nysrc = repmat(srcnorm(2,:),nkappa*nt,1);

  submat = zeros(nkappa,2,nt,2*ns);
  
  submat(:,1,:,1:2:2*ns) = reshape(-coef(1)*(hess(:,:,1).*nxsrc + hess(:,:,2).*nysrc),nkappa,1,nt,ns);
  submat(:,1,:,2:2:2*ns) = reshape(coef(2)*grad(:,:,1),nkappa,1,nt,ns);
    
  submat(:,2,:,1:2:2*ns) = reshape(-coef(1)*(hess(:,:,2).*nxsrc + hess(:,:,3).*nysrc),nkappa,1,nt,ns);
  submat(:,2,:,2:2:2*ns) = reshape(coef(2)*grad(:,:,2),nkappa,1,nt,ns);
  submat = reshape(permute(submat,[1,3,2,4]),nkappa*2*nt,2*ns);
  

otherwise
    error('Unknown quasi-periodic Helmholtz kernel type ''%s''.', type);

end


end