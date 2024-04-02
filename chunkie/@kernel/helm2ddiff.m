function obj = helm2ddiff(type, zks, coefs)
%KERNEL.HELM2DDIFF   Construct the difference of Helmholtz kernels.
%   KERNEL.HELM2DDIFF('s', ZKS) or KERNEL.HELM2DDIFF('single', ZKS) 
%   constructs the difference of single-layer Helmholtz kernels with 
%   wavenumbers ZKS(1) and ZKS(2), scaled by COEFS, i.e.
%   COEFS(1)*KERNEL.HELM2D('s', ZKS(1)) -
%   COEFS(2)*KERNEL.HELM2D('s', ZKS(2))
%
%   KERNEL.HELM2DDIFF('d', ZKS) or KERNEL.HELM2DDIFF('double', ZKS) 
%   constructs the difference of double-layer Helmholtz kernels with 
%   wavenumbers ZKS(1) and ZKS(2), scaled by COEFS, i.e.
%   COEFS(1)*KERNEL.HELM2D('d', ZKS(1)) -
%   COEFS(2)*KERNEL.HELM2D('d', ZKS(2))
%
%   KERNEL.HELM2DDIFF('sp', ZKS) or KERNEL.HELM2DDIFF('sprime', ZKS) 
%   constructs the difference of Neumann data for single-layer Helmholtz
%   kernels with wavenumbers ZKS(1) and ZKS(2), scaled by COEFS, i.e.
%   COEFS(1)*KERNEL.HELM2D('sp', ZKS(1)) -
%   COEFS(2)*KERNEL.HELM2D('sp', ZKS(2))
%
%   KERNEL.HELM2DDIFF('dp', ZKS) or KERNEL.HELM2DDIFF('dprime', ZKS) 
%   constructs the difference of Neumann data for double-layer Helmholtz
%   kernels with wavenumbers ZKS(1) and ZKS(2), scaled by COEFS, i.e.
%   COEFS(1)*KERNEL.HELM2D('dp', ZKS(1)) -
%   COEFS(2)*KERNEL.HELM2D('dp', ZKS(2))
%
%  COEFS is optional. If not given then COEFS is set to [1 1].
%
% See also CHNK.HELM2D.KERN.

if ( nargin < 1 )
    error('Missing Helmholtz kernel type.');
end

if ( nargin < 2 )
    error('Missing Helmholtz wavenumber.');
end

obj = kernel();
obj.name = 'helmholtz';
obj.params.zks = zks;
obj.opdims = [1 1];
if(nargin < 3)
    coefs = [1 1];
end
obj.params.coefs = coefs;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 's') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 's');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 's', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 's', sigma);
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 'd') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 'd');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 'd', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 'd', sigma);
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 'sprime') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 'sprime');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 'sprime', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 'sprime', sigma);
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 'dprime') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 'dprime');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 'dprime', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 'dprime', sigma);
        if ( abs(coefs(1)-coefs(2)) < eps )
            obj.sing = 'log';
        else
            obj.sing = 'hs';
        end

    otherwise
        error('Unknown Helmholtz difference kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end


end
