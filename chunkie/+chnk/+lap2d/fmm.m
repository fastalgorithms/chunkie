function varargout = fmm(eps,srcinfo,targ,type,sigma,varargin)
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
%                type == 'c', combined layer kernel coefs(1)*D + coefs(2)*S
%                type = 'sprime' normal derivative of single layer, 
%                   note that targinfo must contain targ.n
%                type = 'dprime' normal derivative of double layer,
%                   note that targinfo must contain targ.n
%                type = 'stau' tangential derivative of single layer,
%                   note that targinfo must contain targ.d
%   varargin{1} - coef in the combined layer formula, otherwise
%                does nothing
%
% Optional output:
%   pot - potential/neumann data corresponding to the kernel at the target locations
%   grad  - gradient at target locations
%   hess - hessian at target locations
%
%
%

  if isa(targ,'struct')
    targuse = targ.r(:,:);
  else
    targuse = targ;
  end
   pgt = nargout;
   if(nargout == 0)
      warning('LAP2D.FMM: Nothing to compute in LAP2D.FMM, returning with empty array\n');
      return
   end

   srcuse = [];
   mover2pi = -1.0/2/pi;
   srcuse.sources = srcinfo.r(1:2,:);
   if strcmpi(type,'s') || strcmpi(type,'sgrad') || strcmpi(type,'sprime')
      srcuse.charges = mover2pi*sigma(:).';
   end
   if strcmpi(type,'d') || strcmpi(type,'dprime')
      srcuse.dipstr = mover2pi*sigma(:).';
      srcuse.dipvec = srcinfo.n(1:2,:);
   end
   if strcmpi(type,'c')
     coef = varargin{1};
     srcuse.dipstr = coef(1)*mover2pi*sigma(:).';
     srcuse.dipvec = srcinfo.n(1:2,:);
     srcuse.charges = mover2pi*coef(2)*sigma(:).';
   end
   pg = 0;
   pgtuse = min(pgt,3);
   U = rfmm2d(eps,srcuse,pg,targuse,pgtuse);
   if strcmpi(type,'sgrad')
    varargout{1} = U.gradtarg;
    if (pgt >= 2)
     varargout{2} = U.hesstarg([1 2 2 3],:);
    end
    if (pgt >= 3)
        warning('fmm: hessian not available for kernel %s',type)
    end
   else
    varargout{1} = U.pottarg.';
    if strcmpi(type,'sprime') || strcmpi(type,'dprime')
     if isfield(targ,'n')
         varargout{1} = (U.gradtarg(1,:).*targ.n(1,:) + U.gradtarg(2,:).*targ.n(2,:)).'; 
     else
        error('LAP2D.FMM: targets require normal info when evaluating',...
                 'sprime, or dprime'); 
     end
    end

    if strcmpi(type,'stau')
     if isfield(targ,'d')
         ds = sqrt(targ.d(1,:).^2 + targ.d(2,:).^2);
         dx = targ.d(1,:)./ds;
         dy = targ.d(2,:)./ds;
         varargout{1} = (U.gradtarg(1,:).*dx + U.gradtarg(2,:).*dy).'; 
     else
        error('LAP2D.FMM: targets require derivative info when evaluating',...
                 'stau'); 
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
end
