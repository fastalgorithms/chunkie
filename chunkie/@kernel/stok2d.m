function obj = stok2d(type, mu, coefs)
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
%  %
%   KERNEL.STOK2D('cvel', MU) or KERNEL.ELAST2D('cvelocity', MU)
%   constructs the combined-layer Stokes kernel for velocity with viscosity
%   MU. KERNEL.STOK2D('c', MU) and KERNEL.STOK2D('comb', MU) are
%   equivalent.
%
%   KERNEL.STOK2D('cpres', MU) or KERNEL.ELAST2D('cpressure', MU)
%   constructs the combined-layer Stokes kernel for pressure with viscosity
%   MU.
%
%   KERNEL.STOK2D('ctrac', MU) or KERNEL.ELAST2D('ctraction', MU)
%   constructs the combined-layer Stokes kernel for traction with viscosity
%   MU.
%
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
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'svel', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'log';

    case {'spres', 'spressure'}
        obj.type = 'svel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'spres');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'spres', sigma);
        obj.opdims = [1, 2];
	obj.sing = 'pv';

    case {'strac', 'straction'}
        obj.type = 'svel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'strac');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'strac', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'smooth';

    case {'dvel', 'dvelocity', 'd', 'double'}
        obj.type = 'dvel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dvel');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'dvel', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'smooth';

    case {'dpres', 'dpressure'}
        obj.type = 'dpres';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dpres');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'dpres', sigma);
        obj.opdims = [1, 2];
	    obj.sing = 'hs';

    case {'dtrac', 'dtraction'}
        obj.type = 'dtrac';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dtrac');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'dtrac', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'hs';

    case {'cvel', 'cvelocity', 'c', 'comb'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter coefs. Defaulting to 1.');
            coefs = ones(2,1);
        end
        obj.type = 'cvel';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) coefs(1)*chnk.stok2d.kern(mu, s, t, 'dvel') + ...
                          coefs(2)*chnk.stok2d.kern(mu, s, t, 'svel');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'cvel', sigma, coefs);
        obj.opdims = [2, 2];
        obj.sing = 'log';

    case {'cpres', 'cpressure'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter coefs. Defaulting to 1.');
            coefs = ones(2,1);
        end
        obj.type = 'cpres';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) coefs(1)*chnk.stok2d.kern(mu, s, t, 'dpres') + ...
                          coefs(2)*chnk.stok2d.kern(mu, s, t, 'spres');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'cpres', sigma, coefs);
        obj.opdims = [1, 2];
	    obj.sing = 'hs';

    case {'ctrac', 'ctraction'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter coefs. Defaulting to 1.');
            coefs = ones(2,1);
        end
        obj.type = 'ctrac';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) coefs(1)*chnk.stok2d.kern(mu, s, t, 'dtrac') + ...
                          coefs(2)*chnk.stok2d.kern(mu, s, t, 'strac');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'ctrac', sigma, coefs);
        obj.opdims = [2, 2];
        obj.sing = 'hs';

    otherwise
        error('Unknown Stokes kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end


end
