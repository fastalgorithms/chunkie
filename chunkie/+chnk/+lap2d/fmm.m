function varargout = fmm(eps, srcinfo, targinfo, type, sigma, varargin)
%CHNK.LAP2D.FMM   Fast multipole methods for evaluating Laplace layer
%potentials, their gradients, and Hessians.
% 
% Syntax:
%   pot = chnk.lap2d.fmm(eps, srcinfo, targinfo, type, sigma, varargin)
%   [pot, grad] = chnk.lap2d.fmm(eps, srcinfo, targinfo, type, sigma, varargin)
%   [pot, grad, hess] = chnk.lap2d.fmm(eps, srcinfo, targinfo, type, sigma, varargin)
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
%   sigma - density
%   varargin{1} - coef in the combined layer formula, otherwise
%                 does nothing
%
% Optional output:
%   pot - potential/neumann data corresponding to the kernel at the target locations
%   grad  - gradient at target locations
%   hess - hessian at target locations

if ( nargout == 0 )
    warning('CHUNKIE:lap2d:fmm:empty', ...
        'Nothing to compute in LAP2D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r(1:2,:);
switch lower(type)
    case {'s', 'sgrad', 'sprime', 'sp'}
        srcuse.charges = -1/(2*pi)*sigma(:).';
    case {'d', 'dgrad', 'dprime', 'dp'}
        srcuse.dipstr = -1/(2*pi)*sigma(:).';
        srcuse.dipvec = srcinfo.n(1:2,:);
    case {'c', 'cgrad', 'cprime', 'cp'}
        coef = varargin{1};
        srcuse.charges = -1/(2*pi)*coef(2)*sigma(:).';
        srcuse.dipstr  = -1/(2*pi)*coef(1)*sigma(:).';
        srcuse.dipvec  = srcinfo.n(1:2,:);
end

if ( isstruct(targinfo) )
    targuse = targinfo.r(:,:);
else
    targuse = targinfo;
end

pg = 0;
pgt = min(nargout, 3);
switch lower(type)
    case {'sprime', 'dprime', 'cprime','sp','dp','cp','sgrad','dgrad','cgrad','sg','dg','cg'}
        pgt = min(nargout+1,3);
        pgt = max(pgt, 2);
end
U = rfmm2d(eps, srcuse, pg, targuse, pgt);

% Assign potentials
if ( nargout > 0 )
    switch lower(type)
        case {'s', 'd', 'c'}
            varargout{1} = U.pottarg.';
        case {'sgrad', 'dgrad', 'cgrad'}
            varargout{1} = U.gradtarg;
        case {'sprime', 'dprime','cprime','sp','dp','cp'}
            if ( ~isfield(targinfo, 'n') )
                error('CHUNKIE:lap2d:fmm:normals', ...
                    'Targets require normal info when evaluating Laplace kernel ''%s''.', type);
            end
            varargout{1} = ( U.gradtarg(1,:).*targinfo.n(1,:) + ...
                             U.gradtarg(2,:).*targinfo.n(2,:) ).';
        case 'stau'
            if ( ~isfield(targinfo, 'd') )
                error('CHUNKIE:lap2d:fmm:deriv', ...
                    'Targets require derivative info when evaluating Laplace kernel ''%s''.', type);
            end
            ds = sqrt(targinfo.d(1,:).^2 + targinfo.d(2,:).^2);
            dx = targinfo.d(1,:)./ds;
            dy = targinfo.d(2,:)./ds;
            varargout{1} = ( U.gradtarg(1,:).*dx + U.gradtarg(2,:).*dy ).';
        otherwise
            error('CHUNKIE:lap2d:fmm:pot', ...
                'Potentials not supported for Laplace kernel ''%s''.', type);
    end
end

% Assign gradients
if ( nargout > 1 )
    switch lower(type)
        case {'s', 'd', 'c'}
            varargout{2} = U.gradtarg;
        case {'sgrad', 'dgrad', 'cgrad'}
            varargout{2} = U.hesstarg([1 2 2 3],:);
        otherwise
            error('CHUNKIE:lap2d:fmm:grad', ...
                'Gradients not supported for Laplace kernel ''%s''.', type);
    end
end

% Assign Hessians
if ( nargout > 2 )
    switch lower(type)
        case {'s', 'd', 'c'}
            varargout{3} = U.hesstarg;
        otherwise
            error('CHUNKIE:lap2d:fmm:hess', ...
                'Hessians not supported for Laplace kernel ''%s''.', type);
    end
end

end
