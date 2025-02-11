function obj = elast2d(type, lam, mu)
%KERNEL.ELAST2D   Construct the elasticity kernel.
%   KERNEL.ELAST2D('s', LAM, MU) or KERNEL.ELAST2D('single', LAM, MU)
%   constructs the single-layer elasticity kernel with Lamé parameters LAM
%   and MU.
%
%   KERNEL.ELAST2D('strac', LAM, MU) or KERNEL.ELAST2D('straction', LAM,
%   MU) constructs the single-layer elasticity kernel for traction with
%   Lamé parameters LAM and MU.
%
%   KERNEL.ELAST2D('sgrad', LAM, MU) or KERNEL.ELAST2D('straction', LAM,
%   MU) constructs the gradient of the single-layer elasticity kernel with
%   Lamé parameters LAM and MU.
%
%   KERNEL.ELAST2D('d', LAM, MU) or KERNEL.ELAST2D('double', LAM, MU)
%   constructs the double-layer elasticity kernel with Lamé parameters LAM
%   and MU.
%
%   KERNEL.ELAST2D('dalt', LAM, MU) constructs the alternative double-layer
%   elasticity kernel with Lamé parameters LAM and MU.
%
%   KERNEL.ELAST2D('dalttrac', LAM, MU) or 
%   KERNEL.ELAST2D('dalttraction', LAM, MU)
%   constructs the traction of the alternative double-layer elasticity 
%   kernel with Lamé parameters LAM and MU.
%
%   KERNEL.ELAST2D('daltgrad', LAM, MU) constructs the gradient of the 
%   alternative double-layer elasticity kernel with Lamé parameters LAM 
%   and MU.
%
% See also CHNK.ELAST2D.KERN

if ( nargin < 1 )
    error('Missing elasticity kernel type.');
end

if ( nargin < 3 )
    error('Missing Lamé parameters.');
end

obj = kernel();
obj.name = 'elasticity';
obj.params.lam = lam;
obj.params.mu = mu;

obj.opdims = [2 2];

% TODO: Add FMMs.

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.elast2d.kern(lam, mu, s, t, 's');
        obj.sing = 'log';

    case {'sgrad'}
        obj.type = 'sgrad';
        obj.eval = @(s,t) chnk.elast2d.kern(lam, mu, s, t, 'sgrad');
        obj.sing = 'pv';

    case {'strac', 'straction'}
        obj.type = 'strac';
        obj.eval = @(s,t) chnk.elast2d.kern(lam, mu, s, t, 'strac');
        obj.sing = 'pv';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.elast2d.kern(lam, mu, s, t, 'd');
        obj.sing = 'pv';

    case {'dalt'}
        obj.type = 'dalt';
        obj.eval = @(s,t) chnk.elast2d.kern(lam, mu, s, t, 'dalt');
        obj.sing = 'smooth';

    case {'dalttrac','dalttraction'}
        obj.type = 'dalttrac';
        obj.eval = @(s,t) chnk.elast2d.kern(lam, mu, s, t, 'dalttrac');
        obj.sing = 'hs';

    case {'daltgrad'}
        obj.type = 'daltgrad';
        obj.eval = @(s,t) chnk.elast2d.kern(lam, mu, s, t, 'daltgrad');
        obj.sing = 'hs';
        obj.opdims = [4 2];

    otherwise
        error('Unknown elasticity kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end

end
