function varargout = fmm(eps, zk, srcinfo, targinfo, type, sigma, varargin)
%CHNK.HELM2D.FMM   Fast multipole methods for evaluating Helmholtz layer
%potentials, their gradients, and Hessians.
% 
% Syntax:
%   pot = chnk.helm2d.fmm(eps, zk, srcinfo, targinfo, type, sigma, varargin)
%   [pot, grad] = chnk.helm2d.fmm(eps, zk, srcinfo, targinfo, type, sigma, varargin)
%   [pot, grad, hess] = chnk.helm2d.fmm(eps, zk, srcinfo, targinfo, type, sigma, varargin)
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
%   sigma - density
%   varargin{1} - coef in the combined layer formula, otherwise
%                does nothing
%
% Optional Output:
%   pot - potential/neumann data corresponding to the kernel at the target locations
%   grad  - gradient at target locations
%   hess - hessian at target locations

if ( nargout == 0 )
    warning('CHUNKIE:helm2d:fmm:empty', ...
        'Nothing to compute in HELM2D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r(1:2,:);
switch lower(type)
    case {'s', 'sprime'}
        srcuse.charges = sigma(:).';
    case {'d', 'dprime'}
        srcuse.dipstr = sigma(:).';
        srcuse.dipvec = srcinfo.n(1:2,:);
    case 'c'
        coefs = varargin{1};
        srcuse.charges = coefs(2)*sigma(:).';
        srcuse.dipstr  = coefs(1)*sigma(:).';
        srcuse.dipvec  = srcinfo.n(1:2,:);
end

if ( isstruct(targinfo) )
    targuse = targinfo.r(:,:);
else
    targuse = targinfo;
end

pg = 0;
pgt = min(nargout, 2);
U = hfmm2d(eps, zk, srcuse, pg, targuse, pgt);

% Assign potentials
if ( nargout > 0 )
    switch lower(type)
        case {'s', 'd', 'c'}
            varargout{1} = U.pottarg.';
        case {'sprime', 'dprime'}
            if ( ~isfield(targinfo, 'n') )
                error('CHUNKIE:helm2d:fmm:normals', ...
                    'Targets require normal info when evaluating Helmholtz kernel ''%s''.', type);
            end
            varargout{1} = ( U.gradtarg(1,:).*targinfo.n(1,:) + ...
                             U.gradtarg(2,:).*targinfo.n(2,:) ).';
        otherwise
            error('CHUNKIE:lap2d:fmm:pot', ...
                'Potentials not supported for Helmholtz kernel ''%s''.', type);
    end
end

% Assign gradients
if ( nargout > 1 )
    switch lower(type)
        case {'s', 'd', 'c'}
            varargout{2} = U.gradtarg;
        otherwise
            error('CHUNKIE:helm2d:fmm:grad', ...
                'Gradients not supported for Helmholtz kernel ''%s''.', type);
    end
end

% Assign Hessians
if ( nargout > 2 )
    switch lower(type)
        case {'s', 'd', 'c'}
            varargout{3} = U.hesstarg;
        otherwise
            error('CHUNKIE:helm2d:fmm:hess', ...
                'Hessians not supported for Helmholtz kernel ''%s''.', type);
    end
end

end
