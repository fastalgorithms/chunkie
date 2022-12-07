function obj = stok2d(type, mu)
%KERNEL.STOK2D   Construct the Stokes kernel.
%   KERNEL.STOK2D('svel', MU) or KERNEL.ELAST2D('svelocity', MU)
%   constructs the single-layer Stokes kernel for velocity with viscosity
%   MU. KERNEL.STOK2D('s', MU) and KERNEL.STOK2D('single', MU) are
%   equivalent.
%
%   KERNEL.STOK2D('spres', MU) or KERNEL.ELAST2D('spressure', MU)
%   constructs the single-layer Stokes kernel for pressure with viscosity
%   MU.
%
%   KERNEL.STOK2D('strac', MU) or KERNEL.ELAST2D('straction', MU)
%   constructs the single-layer Stokes kernel for traction with viscosity
%   MU.
%
%   KERNEL.STOK2D('dvel', MU) or KERNEL.ELAST2D('dvelocity', MU)
%   constructs the double-layer Stokes kernel for velocity with viscosity
%   MU. KERNEL.STOK2D('d', MU) and KERNEL.STOK2D('double', MU) are
%   equivalent.
%
%   KERNEL.STOK2D('dpres', MU) or KERNEL.ELAST2D('dpressure', MU)
%   constructs the double-layer Stokes kernel for pressure with viscosity
%   MU.
%
%   KERNEL.STOK2D('dtrac', MU) or KERNEL.ELAST2D('dtraction', MU)
%   constructs the double-layer Stokes kernel for traction with viscosity
%   MU.
%
% See also CHNK.STOK2D.KERN.

if ( nargin < 1 )
    error('Missing Stokes kernel type.');
end

if ( nargin < 2 )
    error('Missing Stokes viscosity mu.');
end

obj = kernel();
obj.name = 'stokes';
obj.params.mu = mu;

% TODO: Add FMMs.
% TODO: Add singularities for Stokes.

switch lower(type)

    case {'svel', 'svelocity', 's', 'single'}
        obj.type = 'svel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'svel');
        % mu in denominator

    case {'spres', 'spressure'}
        obj.type = 'svel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'spres');
        % no mu

    case {'strac', 'straction'}
        obj.type = 'svel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'strac');
        % no mu

    case {'dvel', 'dvelocity', 'd', 'double'}
        obj.type = 'dvel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dvel');
        % no mu

    case {'dpres', 'dpressure'}
        obj.type = 'dpres';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dpres');
        % mu in numerator

    case {'dtrac', 'dtraction'}
        obj.type = 'dtrac';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dtrac');
        % mu in numerator

    otherwise
        error('Unknown Stokes kernel type ''%s''.', type);

end

end
