function varargout = fmm(eps,zk,srcinfo,targ,type,sigma,varargin)
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
%                type == 'c', combined layer kernel coefs(1)*D + coefs(2)*S
%                type = 'sprime' normal derivative of single layer, 
%                   note that targinfo must contain targ.n
%                type = 'dprime' normal derivative of double layer,
%                   note that targinfo must contain targ.n
%   varargin{1} - coef in the combined layer formula, otherwise
%                does nothing
%
% Optional Output:
%   pot - potential/neumann data corresponding to the kernel at the target locations
%   grad  - gradient at target locations
%   hess - hessian at target locations
%
%

    if isa(targ,'struct')
        targuse = targ.r(:,:);
    else
        targuse = targ;
    end

   pgt = nargout;
   if(nargout == 0)
      warning('HELM2D.FMM: Nothing to compute in HELM2D.FMM, returning with empty array\n');
      return
   end


   srcuse = [];
   srcuse.sources = srcinfo.r(1:2,:);
   if strcmpi(type,'s') || strcmpi(type,'sprime')
      srcuse.charges = sigma(:).';
   end
   if strcmpi(type,'d') || strcmpi(type, 'dprime')
      srcuse.dipstr = sigma(:).';
      srcuse.dipvec = srcinfo.n(1:2,:);
   end
   if strcmpi(type,'c')
     coefs = varargin{1};
     srcuse.dipstr = coefs(1)*sigma(:).';
     srcuse.dipvec = srcinfo.n(1:2,:);
     srcuse.charges = coefs(2)*sigma(:).';
   end

   pg = 0;

   pgtuse = min(pgt,2);

   U = hfmm2d(eps,zk,srcuse,pg,targuse,pgtuse);
   varargout{1} = U.pottarg.';
   if strcmpi(type,'sprime') || strcmpi(type,'dprime')
     if isfield(targ,'n')
         varargout{1} = (U.gradtarg(1,:).*targ.n(1,:) + U.gradtarg(2,:).*targ.n(2,:)).'; 
     else
        error('HELM2D.FMM: targets require normal info when evaluating',...
                 'sprime, or dprime'); 
     end
   end
   if(pgt>=2) 
       if strcmpi(type, 'sprime') || strcmpi(type,'dprime') || strcmpi(type,'stau')
           warning("Gradients not supported for sprime, dprime, stau")
       else
           varargout{2} = U.gradtarg;
       end
   end
   if(pgt==3)
       if strcmpi(type, 'sprime') || strcmpi(type,'dprime') || strcmpi(type,'stau')
           warning("Hessians not supported for sprime, dprime, stau")
       else
           varargout{3} = U.hesstarg;
       end
   end
end
