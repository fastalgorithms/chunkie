function obj = ostok2d(type, zk, coefs)
%KERNEL.OSTOK2D   Construct the Stokes kernel.
% KERNEL.OSTOK2D('dvel', zk) 
%   constructs the single-layer Oscillatory Stokes kernel for velocity with
%   viscosity zk. KERNEL.OSTOK2D('s', zk) and KERNEL.OSTOK2D('single', zk) 
%   are equivalent.
% KERNEL.OSTOK2D('dvel', zk) 
%   constructs the Double-layer Oscillatory Stokes kernel for velocity with
%   viscosity zk. KERNEL.OSTOK2D('d', zk) and KERNEL.OSTOK2D('double', zk) 
%   are equivalent.
% KERNEL.OSTOK2D('cvel', zk, coefs) 
%   constructs the Combined field Oscillatory Stokes kernel for velocity 
%   with viscosity zk and coefs. KERNEL.OSTOK2D('c', zk, coefs) and 
%   KERNEL.OSTOK2D('comb', zk, coefs) are equivalent. 
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
    case {'dvel', 'dvelocity', 'd', 'double'}
        obj.type = 'dvel';
        obj.eval = @(s,t) chnk.ostok2d.kern(zk, s, t, 'd');
        obj.fmm = [];
        obj.opdims = [2, 2];
        obj.sing = 'log';
    case {'cvel', 'cvelocity', 'c', 'combined'}
        if ( nargin < 2 )
            warning(['Missing combined layer parameter coefs. ' ...
                'Defaulting to [1 1].']);
            coefs = ones(2,1);
        end
        obj.type = 'cvel';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) coefs(1)*chnk.ostok2d.kern(zk, s, t, 'd') + ...
                          coefs(2)*chnk.ostok2d.kern(zk, s, t, 's');
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
