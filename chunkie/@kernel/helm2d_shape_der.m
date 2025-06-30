function obj = helm2d_shape_der(type,zk,vfun)
%KERNEL.helm2D   Construct the helmholtz kernel.
%   KERNEL.helm2D('s') or KERNEL.helm2D('single') constructs the single-layer
%   helmlace kernel.
%
%   KERNEL.helm2D('d') or KERNEL.helm2D('double') constructs the double-layer
%   helmlace kernel.
%
%   KERNEL.helm2D('sp') or KERNEL.helm2D('sprime') constructs the
%   normal derivative of the single-layer helmlace kernel.
%
%   KERNEL.helm2D('sg') or KERNEL.helm2D('sgrad') constructs the
%   gradient of the single-layer helmlace kernel.
%
%   KERNEL.helm2D('dp') or KERNEL.helm2D('dprime') constructs the normal
%   derivative of the double-layer helmlace kernel.
%
%   KERNEL.helm2D('dg') or KERNEL.helm2D('dgrad') constructs the gradient
%   of the double-layer helmlace kernel.
%
%   KERNEL.helm2D('c', coefs) or KERNEL.helm2D('combined', coefs) constructs
%   the combined-layer helmlace kernel with parameter coefs, i.e.,
%   coefs(1)*KERNEL.helm2D('d') + coefs(2)*KERNEL.helm2D('s'). If no
%   value of coefs is specified the default is coefs = [1 1]
%
% See also CHNK.helm2D.KERN.

% author: Dan Fortunato

if ( nargin < 1 )
    error('Missing Helmholtz kernel type.');
end

obj = kernel();
obj.name = 'helmholtz_shape_der';
obj.opdims = [1 1];

switch lower(type)

    case {'dsdz'}
        obj.type = 'dsdz';
        obj.eval = @(s,t) chnk.helm2d_shape_der.kern(s, t, 'dsdz',zk,vfun);
        obj.sing = 'log';

    case {'dspdz'}
        obj.type = 'dspdz';
        obj.eval = @(s,t) chnk.helm2d_shape_der.kern(s, t, 'dspdz',zk,vfun);
        obj.sing = 'log';

    case {'skv'}
        obj.type = 'skv';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'skv',zk,vfun);
        obj.sing = 'log';

    case {'sprimekv'}
        obj.type = 'sprimekv';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'sprimekv',zk,vfun);
        obj.sing = 'smooth';


    case {'vsp', 'vsprime'}
        obj.type = 'vsp';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'vsprime',zk,vfun);
        obj.sing = 'smooth';

    case {'dv'}
        obj.type = 'dv';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'dv',zk,vfun);
        obj.sing = 'smooth';

    case {'vspp-sppv'}
        obj.type = 'vSpp-Sppv';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'vSpp-Sppv',zk,vfun);
        obj.sing = 'pv';

    case {'vdtaustau'}
        obj.type = 'vdtauStau';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'vdtauStau',zk,vfun);
        obj.sing = 'pv';

    case {'vspp_diff'}
        obj.type = 'vSpp-Sppv+vdtauStau';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'vSpp-Sppv+vdtauStau',zk,vfun);
        obj.sing = 'log';

    case {'dpv+sppv'}
        obj.type = 'dpv+sppv';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'dpv+sppv',zk,vfun);
        obj.sing = 'log';


end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end

end


