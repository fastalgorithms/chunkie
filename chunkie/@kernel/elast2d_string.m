function obj = elast2d_string(type, lam, mu,ifimage)
%KERNEL.ELAST2D   Construct the elasticity kernel.
%   KERNEL.ELAST2D('s', LAM, MU) or KERNEL.ELAST2D('single', LAM, MU)
%   constructs the single-layer elasticity kernel with Lamé parameters LAM
%   and MU.
%
%   KERNEL.ELAST2D('strac', LAM, MU) or KERNEL.ELAST2D('straction', LAM,
%   MU) constructs the single-layer elasticity kernel for traction with
%   Lamé parameters LAM and MU.
%
if ( nargin < 1 )
    error('Missing elasticity kernel type.');
end

if ( nargin < 3 )
    error('Missing Lamé parameters.');
end

obj = kernel();
obj.name = 'elasticity_string';
obj.params.lam = lam;
obj.params.mu = mu;

obj.opdims = [2 2];


switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.elast2d_string.kern(lam, mu, s, t, 's',ifimage);
        obj.sing = 'log';

    case {'strac', 'straction'}
        obj.type = 'strac';
        obj.eval = @(s,t) chnk.elast2d_string.kern(lam, mu, s, t, 'strac',ifimage);
        obj.sing = 'removable';

 
    otherwise
        error('Unknown elasticity kernel type ''%s''.', type);

end

    obj.fmm = [];

end
