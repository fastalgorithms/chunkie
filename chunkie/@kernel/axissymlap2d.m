function obj = axissymlap2d(type)
%KERNEL.AXISSYMHELM2D   Construct the axissymmetric Helmholtz kernel.
%   KERNEL.AXISSYMHELM2D('s', ZK) or KERNEL.AXISSYMHELM2D('single', ZK) 
%   constructs the axissymmetric single-layer Helmholtz kernel with 
%   wavenumber ZK.
%
%   KERNEL.AXISSYMHELM2D('d', ZK) or KERNEL.AXISSYMHELM2D('double', ZK) 
%   constructs the axissymmetric double-layer Helmholtz kernel with 
%   wavenumber ZK.
%
%   KERNEL.AXISSYMHELM2D('sp', ZK) or KERNEL.AXISSYMHELM2D('sprime', ZK) 
%   constructs the normal derivative of the axissymmetric single-layer 
%   Helmholtz kernel with wavenumber ZK.
%
%   KERNEL.AXISSYMHELM2D('c', ZK, COEFS) or 
%   KERNEL.AXISSYMHELM2D('combined', ZK, COEFS)
%   constructs the combined-layer axissymmetric Helmholtz kernel with 
%   wavenumber ZK and parameter COEFS, 
%   i.e., COEFS(1)*KERNEL.AXISSYMHELM2D('d', ZK) + 
%         COEFS(2)*KERNEL.AXISSYMHELM2D('s', ZK).
%
% NOTES: The axissymetric kernels are currently supported only for purely
% real or purely imaginary wave numbers
%
% See also CHNK.AXISSYMHELM2D.KERN.

if ( nargin < 1 )
    error('Missing Laplace kernel type.');
end

obj = kernel();
obj.name = 'axissymlaplace';
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 's');
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 's');
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'd');
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'd');
        obj.fmm = [];
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'sprime');
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'sprime');
        obj.fmm = [];
        obj.sing = 'log';

    otherwise
        error('Unknown axissym Laplace kernel type ''%s''.', type);

end

end
