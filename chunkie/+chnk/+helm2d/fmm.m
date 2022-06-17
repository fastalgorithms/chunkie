function [pot,varargout] = fmm(eps,zk,srcinfo,targ,type,sigma,pgt,varargin)
%CHNK.HELM2D.FMM fast multipole methods for evaluating helmholtz layer
%potentials, their gradients, and hessians
% 
% Syntax: pot = chnk.helm2d.fmm(zk,srcinfo,targinfo,sigma,type,pg,varargin)
% [pot,grad] = chnk.helm2d.fmm(zk,srcinfo,targinfo,sigma,type,pg,varargin)
% [pot,grad,hess] = chnk.helm2d.fmm(zk,srcinfo,targinfo,sigma,type,pg,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). 
%  
% Kernels based on G(x,y) = i/4 H_0^{(1)}(zk |x-y|)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
%
% Note: the density should be scaled by the weights
%
% Input:
%   eps - precision requested
%   zk - complex number, Helmholtz wave number
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
%                type == 'c', combined layer kernel D + i eta S
%   varargin{1} - eta in the combined layer formula, otherwise
%                does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2D.FMM
%
   srcuse = [];
   srcuse.sources = srcinfo.r(1:2,:);
   if strcmpi(type,'s')
      srcuse.charges = sigma(:).';
   end
   if strcmpi(type,'d')
      srcuse.dipstr = sigma(:).';
      srcuse.dipvec = srcinfo.n(1:2,:);
   end
   if strcmpi(type,'c')
     eta = varargin{1};
     srcuse.dipstr = sigma(:).';
     srcuse.dipvec = srcinfo.n(1:2,:);
     srcuse.charges = eta*sigma(:).';
   end
   pg = 0;
   U = hfmm2d(eps,zk,srcuse,pg,targ,pgt);
   pot = U.pottarg.';
   if(pgt>=2) 
       varargout{1} = U.gradtarg; 
   end
   if(pgt==3) 
       varargout{2} = U.hesstarg; 
   end
end
