function obj = ostok2d(type, zk)
%KERNEL.STOK2D   Construct the Stokes kernel.
%   KERNEL.STOK2D('svel', MU) or KERNEL.ELAST2D('svelocity', MU)
%   constructs the single-layer Stokes kernel for velocity with viscosity
%   MU. KERNEL.STOK2D('s', MU) and KERNEL.STOK2D('single', MU) are
%   equivalent.
%
%
% See also CHNK.OSTOK2D.KERN.

% author: Kshitij Sinha

if ( nargin < 1 )
    error('Missing Oscillatory Stokes kernel type.');
end

if ( nargin < 2 )
    error('Missing Oscillatory Stokes wave number k.');
end

obj = kernel();
obj.name = 'ostokes';
obj.params.zk = zk;


switch lower(type)

    case {'svel', 'svelocity', 's', 'single'}
        obj.type = 'svel';
        obj.eval = @(s,t) chnk.ostok2d.kern(zk, s, t, 's');
        obj.fmm = [];
        obj.opdims = [2, 2];
        obj.sing = 'log';
    otherwise
        error('Unknown Oscillatory Stokes kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end

end
