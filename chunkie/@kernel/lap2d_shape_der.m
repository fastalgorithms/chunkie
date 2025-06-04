function obj = lap2d_shape_der(type,vfun)
%KERNEL.LAP2D   Construct the Laplace kernel.
%   KERNEL.LAP2D('s') or KERNEL.LAP2D('single') constructs the single-layer
%   Laplace kernel.
%
%   KERNEL.LAP2D('d') or KERNEL.LAP2D('double') constructs the double-layer
%   Laplace kernel.
%
%   KERNEL.LAP2D('sp') or KERNEL.LAP2D('sprime') constructs the
%   normal derivative of the single-layer Laplace kernel.
%
%   KERNEL.LAP2D('sg') or KERNEL.LAP2D('sgrad') constructs the
%   gradient of the single-layer Laplace kernel.
%
%   KERNEL.LAP2D('dp') or KERNEL.LAP2D('dprime') constructs the normal
%   derivative of the double-layer Laplace kernel.
%
%   KERNEL.LAP2D('dg') or KERNEL.LAP2D('dgrad') constructs the gradient
%   of the double-layer Laplace kernel.
%
%   KERNEL.LAP2D('c', coefs) or KERNEL.LAP2D('combined', coefs) constructs
%   the combined-layer Laplace kernel with parameter coefs, i.e.,
%   coefs(1)*KERNEL.LAP2D('d') + coefs(2)*KERNEL.LAP2D('s'). If no
%   value of coefs is specified the default is coefs = [1 1]
%
% See also CHNK.LAP2D.KERN.

% author: Dan Fortunato

if ( nargin < 1 )
    error('Missing Laplace kernel type.');
end

obj = kernel();
obj.name = 'laplace_shape_der';
obj.opdims = [1 1];

switch lower(type)

    case {'dsdz'}
        obj.type = 'dsdz';
        obj.eval = @(s,t) chnk.lap2d_shape_der.kern(s, t, 'dsdz',vfun);
        obj.sing = 'log';

    case {'dspdz'}
        obj.type = 'dspdz';
        obj.eval = @(s,t) chnk.lap2d_shape_der.kern(s, t, 'dspdz',vfun);
        obj.sing = 'log';

    case {'skv'}
        obj.type = 'skv';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'skv',vfun);
        obj.sing = 'log';

    case {'sprimekv'}
        obj.type = 'sprimekv';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'sprimekv',vfun);
        obj.sing = 'smooth';


    case {'vsp', 'vsprime'}
        obj.type = 'vsp';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'vsprime',vfun);
        obj.sing = 'smooth';

    case {'dv'}
        obj.type = 'dv';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'dv',vfun);
        obj.sing = 'smooth';

    case {'vspp-sppv'}
        obj.type = 'vSpp-Sppv';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'vSpp-Sppv',vfun);
        obj.sing = 'pv';

    case {'vdtaustau'}
        obj.type = 'vdtauStau';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'vdtauStau',vfun);
        obj.sing = 'pv';

    case {'vspp_diff'}
        obj.type = 'vSpp-Sppv+vdtauStau';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'vSpp-Sppv+vdtauStau',vfun);
        obj.sing = 'log';

    case {'dpv+sppv'}
        obj.type = 'dpv+sppv';
        obj.eval = @(s,t) chnk.lap2dv.kern(s, t, 'dpv+sppv',vfun);
        obj.sing = 'log';


end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end

end


