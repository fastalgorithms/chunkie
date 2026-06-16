function submat = kern(srcinfo,targinfo,type,kappa,d,s0,sn,l,ising,nsub,varargin)
%CHNK.LAP2DQUAS.KERN quasiperiodic Laplace layer potential kernels in 2D
%
% Syntax: submat = chnk.lap2dquas.kern(srcinfo,targinfo,type,kappa,d,s0,sn,l,ising)
%         submat = chnk.lap2dquas.kern(srcinfo,targinfo,type,kappa,d,s0,sn,l,ising,nsub)
%
% Let x be targets and y be sources, with n_x and n_y the corresponding
% unit outward normals and tau_x, tau_y the unit tangents.
%
% Kernel definitions (based on the quasiperiodic Laplace Green's function G):
%   S(x,y)  = G(x,y)
%   D(x,y)  = grad_{n_y} G(x,y)
%   S'(x,y) = grad_{n_x} G(x,y)
%   D'(x,y) = grad_{n_x} grad_{n_y} G(x,y)
%
% Input:
%   srcinfo  - source descriptor in ptinfo struct format:
%                srcinfo.r  - positions (2,:)
%                srcinfo.n  - unit outward normal (2,:) [needed for 'd','hilb','hilbprime','dp']
%   targinfo - target descriptor in ptinfo struct format:
%                targinfo.r  - positions (2,:)
%                targinfo.n  - unit outward normal (2,:) [needed for 'sp','stau','hilbprime','dp']
%   type     - string, determines kernel type:
%                's'  or 'single'    - single layer S
%                'd'  or 'double'    - double layer D
%                'sp' or 'sprime'    - normal derivative of single layer S'
%                'st' or 'stau'      - tangential derivative of single layer
%                'hilb'              - quasiperiodic Hilbert transform
%                                      (2 * tangential adjoint of S)
%                'hilbprime'         - tangential derivative of Hilbert transform
%                'dp' or 'dprime'    - normal derivative of double layer D'
%   kappa    - (nkappa,1) array of quasiperiodic phase parameters
%   d        - period (scalar)
%   s0       - (nkappa,1) n=0 lattice sum constant
%                  (see chnk.lap2dquas.latticecoefs)
%   sn       - (nkappa, N) lattice sum coefficients for orders 1..N
%                  (see chnk.lap2dquas.latticecoefs)
%   l        - number of explicit periodic copies on each side
%   ising    - if 1, include the free-space singular part; if 0, periodic
%                  images only (smooth kernel)
%   nsub     - (optional, default 0) number of additional source copies to
%                  subtract near the source
%
% Output:
%   submat - (nkappa*ntarg, nsrc) kernel matrix; rows correspond to targets,
%            columns to sources
%
% see also CHNK.LAP2DQUAS.GREEN, CHNK.LAP2DQUAS.LATTICECOEFS

src = srcinfo.r(:,:);
targ = targinfo.r(:,:);

[~,ns] = size(src);
[~,nt] = size(targ);
nkappa = length(kappa);

if nargin < 10
    nsub = 0;
end

switch lower(type)
% double layer
case {'d', 'double'}
    %srcnorm = chnk.normal2d(srcinfo);
    srcnorm = srcinfo.n(:,:);
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising,nsub);
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

% normal derivative of single layer
case {'sp', 'sprime'}
    targnorm = targinfo.n;
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    nx = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nx = reshape(nx,[],ns);
    ny = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    ny = reshape(ny,[],ns);
    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);


% Tangential derivative of single layer
case {'stau'}
    targnorm = targinfo.n;
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    % nx = repmat((targnorm(1,:)).',1,ns);
    % ny = repmat((targnorm(2,:)).',1,ns);
    nx = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nx = reshape(nx,[],ns);
    ny = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    ny = reshape(ny,[],ns);

    submat = (-grad(:,:,1).*ny + grad(:,:,2).*nx);

% Hilbert transform (two times the adjoint of stau)
case {'hilb'} 
    srcnorm = srcinfo.n;
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising,nsub);
    % nx = repmat((srcnorm(1,:)),nt,1);
    % ny = repmat((srcnorm(2,:)),nt,1);
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);

    submat = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

% tangential derivative Hilbert transform (tau of two times the adjoint of stau)
case {'hilbprime'} 
    srcnorm = srcinfo.n;
    targnorm = targinfo.n;
    [~,~,hess] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising,nsub);
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);


    submat = 2*(-(hess(:,:,1).*ny - hess(:,:,2).*nx).*nytarg ...
        + (hess(:,:,2).*ny - hess(:,:,3).*nx).*nxtarg);

% normal derivative of double layer
case {'dprime','dp'}
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    nxsrc = repmat(srcnorm(1,:),nkappa*nt,1);
    nysrc = repmat(srcnorm(2,:),nkappa*nt,1);
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);
    submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);

% single layer
case {'s', 'single'}
    submat = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,1);

otherwise
    error('Unknown quasiperiodic Laplace kernel type ''%s''.', type);
end

end

