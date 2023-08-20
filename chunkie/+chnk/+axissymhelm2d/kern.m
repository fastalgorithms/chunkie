function submat = kern(zk, srcinfo, targinfo, type, varargin)
%CHNK.AXISSYMHELM2D.KERN axissymmetric Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnk.axissymhelm2d.kern(zk,srcinfo,targingo,type,htables)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%
% Here the first and second components correspond to the r and z
% coordinates respectively. 
%
% Kernels based on G(x,y) = \int_{0}^{\pi} e^{i k d(t)}/(d(t)) \, dt \, 
% where d(t) = \sqrt(r^2 + r'^2 - 2rr' \cos(t) + (z-z')^2) with
% x = (r,z), and y = (r',z')
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number, must be purely real, or purely 
%        imaginary, doesn't support other complex wavenumbers yet
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
%                type == 'ddiff', D_{k}' - D_{ik}', for this routine k must be real
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
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
%

src = srcinfo.r; 
targ = targinfo.r;

[~, ns] = size(src);
[~, nt] = size(targ);

if strcmpi(type, 'd')
    srcnorm = srcinfo.n;
    [~, grad] = chnk.axissymhelm2d.green(zk, src, targ);
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);
    % Due to lack of translation invariance in r, no sign flip needed, 
    % as gradient is computed with repsect to r'
    submat = (grad(:,:,2).*nx - grad(:,:,3).*ny);
end

if strcmpi(type, 'sprime')
    targnorm = targinfo.n;
    [~, grad] = chnk.axissymhelm2d.green(zk, src, targ);

    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    submat = (grad(:,:,1).*nx + grad(:,:,3).*ny);
end

if strcmpi(type, 's')
    submat = chnk.axissymhelm2d.green(zk, src, targ);
end


if strcmpi(type, 'sdiff')
    ifdiff = 1;
    submat = chnk.axissymhelm2d.green(zk, src, targ, ifdiff);
end

if strcmpi(type, 'c')
    srcnorm = srcinfo.n; 
    coef = ones(2,1);
    if (nargin == 5); coef = varargin{1}; end
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);
    [submats, grad] = chnk.axissymhelm2d.green(zk, src, targ);
    % Due to lack of translation invariance in r, no sign flip needed, 
    % as gradient is computed with respect to r'
    submat = coef(1)*(grad(:,:,2).*nx - grad(:,:,3).*ny) + coef(2)*submats;
end


if strcmpi(type, 'dprime')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  [~,~,hess] = chnk.axissymhelm2d.green(zk,src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ... 
              - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;
end



if strcmpi(type, 'dprimediff')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  ifdiff = 1;
  [~,~,hess] = chnk.axissymhelm2d.green(zk,src,targ,ifdiff);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ...
      - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;
end

