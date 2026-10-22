function varargout = fmm(eps, mu, srcinfo, targinfo, type, sigma, varargin)
%CHNK.stok2D.FMM   Fast multipole methods for evaluating Stokes layer
%potentials, their gradients, and Hessians.
% 
% Syntax:
%   [vel, pres, grad, presgrad] = 
%     chnk.stok2d.fmm(eps, mu, srcinfo, targinfo, '*vel', sigma, varargin)
%   [pres, presgrad] = 
%     chnk.stok2d.fmm(eps, mu, srcinfo, targinfo, '*pres', sigma, varargin)
%   trac = 
%     chnk.stok2d.fmm(eps, mu, srcinfo, targinfo, '*trac', sigma, varargin)
%
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined).
%  
% Kernels based on G(x,y) = 
%    1/4\pi \mu( r_{i} r_{j}/r^2 - \log{|x-y|} \delta_{ij}) 
% 
% Note: the density should be scaled by the weights
% Note: pressure gradients are currently unavailable
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
%                type == 'dvel', double layer kernel D
%                type == 'svel', single layer kernel S
%                type == 'cvel', combined layer kernel coefs(1)*D + coefs(2)*S
%                type = 'spres', pressure corresponding to single layer
%                type = 'dpres', pressure corresponding to double layer
%                type = 'cpres', pressure corresponding to combined layer
%                type = 'strac' traction corresponding to single layer,
%                   note that targinfo must contain targ.n
%                type = 'dtrac' traction corresponding to double layer,
%                   note that targinfo must contain targ.n
%                type = 'ctrac' traction corresponding to combined layer,
%                   note that targinfo must contain targ.n
%                type = 'sgrad' gradient corresponding to single layer,
%                type = 'dgrad' gradient corresponding to double layer,
%                type = 'cgad' gradient corresponding to combined layer,
%   sigma - density, (2*ns,ndens) for ndens densities
%   varargin{1} - coef in the combined layer formula, otherwise
%                 does nothing
%
% Optional output:
%   vel  - velocity corresponding to kernel at target locations
%   pres - pressure at target locations
%   trac - traction at target locations

if ( nargout == 0 )
    warning('CHUNKIE:stok2d:fmm:empty', ...
        'Nothing to compute in stok2D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r(1:2,:);
[m, ns] = size(srcuse.sources);
ndens = numel(sigma)/(m*ns);
% the densities as (ndens,2,ns), as stfmm2d expects for nd = ndens
sigma = permute(reshape(sigma, [m, ns, ndens]), [3 1 2]);
srcuse.nd = ndens;
if ~any(strcmpi(type, {'svel', 'strac', 'spres', 'sgrad'}))
    strsvec = repmat(reshape(srcinfo.n(1:2,:), [1 m ns]), [ndens 1 1]);
end
switch lower(type)
    case {'svel', 'strac', 'spres', 'sgrad'}
        srcuse.stoklet = 1/(2*pi*mu)*sigma;
    case {'dvel', 'dtrac', 'dpres', 'dgrad'}
        srcuse.strslet = -1/(2*pi)*sigma;
        srcuse.strsvec = strsvec;
    case {'cvel', 'ctrac', 'cpres', 'cgrad'}
        coef = varargin{1};
        srcuse.stoklet = 1/(2*pi*mu)*coef(2)*sigma;
        srcuse.strslet  = -1/(2*pi)*coef(1)*sigma;
        srcuse.strsvec  = strsvec;
end

if ( isstruct(targinfo) )
    targuse = targinfo.r(:,:);
else
    targuse = targinfo;
end

ifppreg = 0;
ifppregtarg = 1;
switch lower(type)
    case {'svel', 'dvel', 'cvel'}
        ifppregtarg = 1;
    case {'spres', 'dpres', 'cpres'}
        ifppregtarg = 2;
    case {'strac', 'dtrac', 'ctrac', 'sgrad', 'dgrad', 'cgrad'}
        ifppregtarg = 3;
end

if nargout > 1
    ifppregtarg = 2;
end

if nargout > 2
    ifppregtarg = 3;
end

U = stfmm2d(eps, srcuse, ifppreg, targuse, ifppregtarg);
[m, nt] = size(targuse);
% the outputs with the densities last: velocity (2,nt,ndens), pressure
% (nt,ndens) and gradient (2,2,nt,ndens)
pot = permute(reshape(U.pottarg, [ndens m nt]), [2 3 1]);
if ( ifppregtarg >= 2 )
    pre = reshape(U.pretarg, [ndens nt]).';
end
if ( ifppregtarg >= 3 )
    grad = permute(reshape(U.gradtarg, [ndens m m nt]), [2 3 4 1]);
end
% Assign potentials
if ( nargout > 0 )
    switch lower(type)
        case {'svel', 'dvel', 'cvel'}
            varargout{1} = reshape(pot, [m*nt, ndens]);
            if nargout > 1
                varargout{2} = pre*mu;
            end
            if nargout > 2
                varargout{3} = grad;
            end
            if nargout > 3
                str_err = ['CHUNKIE:stok2d:fmm: too many output ',...
                    'arguments for Stokes kernel ''%s''.'];
                error(str_err, type);
            end
        case {'spres', 'dpres', 'cpres'}
            varargout{1} = pre*mu;
            if nargout > 1
                str_err = ['CHUNKIE:stok2d:fmm: too many output ',...
                    'arguments for Stokes kernel ''%s''.'];
                error(str_err, type);
            end
        case {'sgrad', 'dgrad', 'cgrad'}
            varargout{1} = reshape(grad, [m*m*nt, ndens]);
            if nargout > 1
                str_err = ['CHUNKIE:stok2d:fmm: too many output ',...
                    'arguments for Stokes kernel ''%s''.'];
                error(str_err, type);
            end

        case {'strac', 'dtrac', 'ctrac'}
            if ( ~isfield(targinfo, 'n') )
                error('CHUNKIE:stok2d:fmm:normals', ...
                    'Targets require normal info when evaluating Stokes kernel ''%s''.', type);
            end
            eu = grad + permute(grad, [2 1 3 4]);
            euxx = reshape(eu(1,1,:,:), [nt ndens]);
            euxy = reshape(eu(1,2,:,:), [nt ndens]);
            euyy = reshape(eu(2,2,:,:), [nt ndens]);
            f = zeros(m*nt,ndens);
            p = pre;
            ntx = targinfo.n(1,:).';
            nty = targinfo.n(2,:).';
            f(1:2:end,:) = -p.*ntx + euxx.*ntx + euxy.*nty;
            f(2:2:end,:) = -p.*nty + euxy.*ntx + euyy.*nty;
            varargout{1} = mu*f;
            if nargout > 1
                str_err = ['CHUNKIE:stok2d:fmm: too many output ',...
                    'arguments for Stokes kernel ''%s''.'];
                error(str_err, type);
            end
        otherwise
            error('CHUNKIE:stok2d:fmm:pot', ...
                'Potentials not supported for Stokes kernel ''%s''.', type);
    end
end



end
