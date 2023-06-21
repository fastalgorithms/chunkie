function [pot,varargout] = fmm(eps,srcinfo,targ,type,sigma,pgt,varargin)
%CHNK.LAP2D.FMM fast multipole methods for evaluating Laplace layer
%potentials, their gradients, and hessians
% 
% Syntax: pot = chnk.lap2d.fmm(srcinfo,targinfo,sigma,type,pg,varargin)
% [pot,grad] = chnk.lap2d.fmm(srcinfo,targinfo,sigma,type,pg,varargin)
% [pot,grad,hess] = chnk.helm2d.fmm(srcinfo,targinfo,sigma,type,pg,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). 
%  
% Kernels based on G(x,y) = -1/2\pi \log{|x-y|}
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
%
% Note: the density should be scaled by the weights
%
% Input:
%   eps - precision requested
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%                ptinfo.n - normal (2,:) array
%                     
%   targinfo - description of targets in ptinfo struct format,
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'c', combined layer kernel D + eta S
%   varargin{1} - eta in the combined layer formula, otherwise
%                does nothing
%
% Output:
%   pot - potential corresponding to the kernel at the target locations
%
% Optional output
%   grad  - gradient at target locations
%   hess - hessian at target locations
%
%

  if isa(targ,'struct')
    targuse = targ.r(:,:);
  else
    targuse = targ;
  end
   srcuse = [];
   mover2pi = -1.0/2/pi;
   srcuse.sources = srcinfo.r(1:2,:);
   if strcmpi(type,'s')
      srcuse.charges = mover2pi*sigma(:).';
   end
   if strcmpi(type,'sgrad')
      srcuse.charges = mover2pi*sigma(:).';
   end
   if strcmpi(type,'d')
      srcuse.dipstr = mover2pi*sigma(:).';
      srcuse.dipvec = srcinfo.n(1:2,:);
   end
   if strcmpi(type,'c')
     eta = varargin{1};
     srcuse.dipstr = mover2pi*sigma(:).';
     srcuse.dipvec = srcinfo.n(1:2,:);
     srcuse.charges = mover2pi*eta*sigma(:).';
   end
   pg = 0;
   U = rfmm2d(eps,srcuse,pg,targuse,pgt);
   if strcmpi(type,'sgrad')
    pot = U.gradtarg;
    if (pgt >= 2)
     varargout{1} = U.hesstarg([1 2 2 3],:);
    end
    if (pgt >= 3)
        warning('fmm: hessian not available for kernel %s',type)
    end
   else
    pot = U.pottarg.';
    if(pgt>=2) 
       varargout{1} = U.gradtarg; 
    end
    if(pgt==3) 
       varargout{2} = U.hesstarg; 
    end
   end
end
