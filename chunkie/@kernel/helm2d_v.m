function obj = helm2d_v(type,zk,vfun)

if ( nargin < 1 )
    error('Missing Helmholtz shape kernel type.');
end

obj = kernel();
obj.name = 'helmholtz_shape';
obj.opdims = [1 1];

switch lower(type)

    case {'spdk'}
        obj.type = 'spdk';
        obj.eval = @(s,t) chnk.helm2dv.kern(s, t, 'spdk',zk);
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


