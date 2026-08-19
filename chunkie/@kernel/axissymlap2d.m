function obj = axissymlap2d(type, n)
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
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 's', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 's', n);
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'd', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'd', n);
        obj.fmm = [];
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'sprime', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'sprime', n);
        obj.fmm = [];
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'dprime', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'dprime', n);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'spp', 'sprimeprime'}
        obj.type = 'spp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'sprimeprime', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'sprimeprime', n);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'q'}
        obj.type = 'q';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'q', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'q', n);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'q_sum_dp'}
        obj.type = 'q_sum_dp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'q_sum_dp', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'q_sum_dp', n);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'spp_sum_dp'}
        obj.type = 'spp_sum_dp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern(s, t, [0,0], 'spp_sum_dp', n);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern(s, t, o, 'spp_sum_dp', n);
        obj.fmm = [];
        obj.sing = 'hs';
     
    otherwise
        error('Unknown axissym Laplace kernel type ''%s''.', type);

end

end
