function obj = axissymlap2d(type, m)
%KERNEL.AXISSYMLAP2D   Construct the axissymmetric Laplace kernel.
%   KERNEL.AXISSYMLAP2D('s', m) or KERNEL.AXISSYMLAP2D('single', m) 
%   constructs the axissymmetric single-layer Laplace kernel for the
%   mth Fourier mode (m >= 0)
%
%   KERNEL.AXISSYMLAP2D('d', m) or KERNEL.AXISSYMLAP2D('double', m) 
%   constructs the axissymmetric double-layer Laplace kernel for the 
%   mth Fourier mode (m >= 0)
%
%   KERNEL.AXISSYMLAP2D('sp', m) or KERNEL.AXISSYMLAP2D('sprime', m) 
%   constructs the normal derivative of the axissymmetric single-layer 
%   Laplace kernel for the mth Fourier mode (m >= 0)
%
%   KERNEL.AXISSYMLAP2D('dp', m) or KERNEL.AXISSYMLAP2D('dprime', m) 
%   constructs the normal derivative of the axissymmetric double-layer 
%   Laplace kernel for the mth Fourier mode (m >= 0)
%
% if m is not provided, it defaults to m=0.
%
% See also CHNK.AXISSYMHELM2D.KERN.

if ( nargin < 1 )
    error('Missing Laplace kernel type.');
end
if ( nargin < 2 )
    m = 0;
end

obj = kernel();
obj.name = 'axissymlaplace';
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 's',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 's',m);
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 'd',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 'd',m);
        obj.fmm = [];
        obj.sing = 'smooth';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 'sprime',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 'sprime',m);
        obj.fmm = [];
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 'dprime',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 'dprime',m);
        obj.fmm = [];
        obj.sing = 'hs';

    otherwise
        error('Unknown axissym Laplace kernel type ''%s''.', type);

end

end
