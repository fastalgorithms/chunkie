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
%   KERNEL.HELM2DDIFF('c', ZKS) or KERNEL.HELM2DDIFF('combined', ZKS) 
%   constructs the difference of combined-layer Helmholtz kernels with 
%   wavenumbers ZKS(1) and ZKS(2), scaled by COEFS, i.e.
%   KERNEL.HELM2D('c', ZKS(1), COEFS(:,1)) -
%   KERNEL.HELM2D('c', ZKS(2), COEFS(:,2))
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
%   KERNEL.HELM2DDIFF('cp', ZKS) or KERNEL.HELM2DDIFF('cprime', ZKS) 
%   constructs the difference of Neumann for combined-layer Helmholtz 
%   kernels with wavenumbers ZKS(1) and ZKS(2), scaled by COEFS, i.e.
%   KERNEL.HELM2D('cp', ZKS(1), COEFS(:,1)) -
%   KERNEL.HELM2D('cp', ZKS(2), COEFS(:,2))
%
%   KERNEL.HELM2DDIFF('all', ZKS, coefs) or 
%   KERNEL.HELM2DDIFF('tsys', ZKS, coefs)
%   constructs the (2x2) matrix of difference kernels of
%   [D, S; D' S'] scaled by coefs where coefs is (2x2x2) tensor, i.e.
%   D part is scaled as COEFS(1,1,1)*KERNEL.HELM2D('d', zks(1)) - 
%   COEFS(1,1,2)*KERNEL.HELM2D('d', ZKS(2)), and so on.
%
%   KERNEL.HELM2DDIFF('trans_rep', ZK, COEFS) or KERNEL.HELM2DDIFF('trep', ZK, COEFS) 
%   constructs the difference of transmission repretsentations, i.e. the (1x2) matrix of
%   kernels [D, S] scaled by coefs where coefs is (1x2) matrix, i.e. D part
%   is scaled as COEFS(1,1)*KERNEL.HELM2DDIFF('d', ZKS(1)) - 
%   COEFS(1,2)*KERNEL.HELM2DDIFF('d', ZKS(2)), and so on.
%
%   KERNEL.HELM2DDIFF('trans_rep_p', ZK, COEFS) or KERNEL.HELM2DDIFF('trep_p', ZK, COEFS) 
%   constructs the difference of derivative of the transmission repretsentation, i.e. the
%   (1x2) matrix of kernels [D', S'] scaled by coefs where coefs is (1x2)
%   matrix, i.e. D part is scaled as COEFS(1,1)*KERNEL.HELM2DDIFF('dp', ZKS(1)) - 
%   COEFS(1,2)*KERNEL.HELM2DDIFF('dp', ZKS(2)), and so on.
%
%   KERNEL.HELM2DDIFF('c2trans', ZK, COEFS) or KERNEL.HELM2DDIFF('c2t', ZK, COEFS) 
%   evaluates the difference of combined-layer Helmholtz kernel and its derivative, i.e.
%   the (2x1) matrix of kernels [C; C'] scaled by coefs where coefs is
%   (2x2) matrix, i.e. kernel returns 
%   KERNEL.HELM2DDIFF('c2trans',zks(1),coefs(:,:,1)) - 
%   KERNEL.HELM2DDIFF('c2trans',zks(2),coefs(:,:,2)) - 
%
%  COEFS is optional. If not given then COEFS is set to all ones.
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
    if strcmpi(type, 'all') || strcmpi(type, 'trans_sys') || strcmpi(type, 'tsys')
      coefs = ones(2,2,2);
    end
end
obj.params.coefs = coefs;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 's_diff') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 's_diff');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 's', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 's', sigma);
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 'd_diff') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 'd_diff');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 'd', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 'd', sigma);
        obj.sing = 'log';

    case {'c', 'combined'}
        if (nargin < 3)
            coefs = ones(2,2);
        end
        obj.type = 'c';
        obj.eval = @(s,t) chnk.helm2d.kern(zks(1), s, t, 'c_diff',coefs(:,1)) - ...
                          chnk.helm2d.kern(zks(2), s, t, 'c_diff',coefs(:,2));
        obj.fmm  = @(eps,s,t,sigma) ...
            chnk.helm2d.fmm(eps, zks(1), s, t, 'c', sigma,coefs(:,1)) - ...
            chnk.helm2d.fmm(eps, zks(2), s, t, 'c', sigma,coefs(:,2));
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 'sprime_diff') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 'sprime_diff');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 'sprime', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 'sprime', sigma);
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) coefs(1)*chnk.helm2d.kern(zks(1), s, t, 'dprime_diff') - ...
                          coefs(2)*chnk.helm2d.kern(zks(2), s, t, 'dprime_diff');
        obj.fmm  = @(eps,s,t,sigma) ...
            coefs(1)*chnk.helm2d.fmm(eps, zks(1), s, t, 'dprime', sigma) - ...
            coefs(2)*chnk.helm2d.fmm(eps, zks(2), s, t, 'dprime', sigma);
        if ( abs(coefs(1)-coefs(2)) < eps )
            obj.sing = 'log';
        else
            obj.sing = 'hs';
        end

    case {'cp', 'cprime'}
        if (nargin < 3)
            coefs = ones(2,2);
        end
        obj.type = 'c';
        obj.eval = @(s,t) chnk.helm2d.kern(zks(1), s, t, 'cp_diff',coefs(:,1)) - ...
                          chnk.helm2d.kern(zks(2), s, t, 'cp_diff',coefs(:,2));
        obj.fmm  = @(eps,s,t,sigma) ...
            chnk.helm2d.fmm(eps, zks(1), s, t, 'cp', sigma,coefs(:,1)) - ...
            chnk.helm2d.fmm(eps, zks(2), s, t, 'cp', sigma,coefs(:,2));
        obj.sing = 'log';

    case {'all', 'trans_sys', 'tsys'}
        if (nargin < 3)
            coefs = ones(2,2,2);
        end
        obj.type = 'all';
        obj.opdims = [2,2];
        obj.eval = @(s,t) chnk.helm2d.kern(zks(1), s, t, 'all_diff', coefs(:,:,1)) - ...
                          chnk.helm2d.kern(zks(2), s, t, 'all_diff', coefs(:,:,2));
        obj.fmm  = []; 
        if ( abs(coefs(2,1,1)-coefs(2,1,2)) < eps )
            obj.sing = 'log';
        else
            obj.sing = 'hs'; 
        end
        
    case {'trans_rep','trep'} 
        if ( nargin < 3 )
            coefs = ones(2,2);
        end
        obj.type = 'trep';
        obj.opdims = [1,2];
        obj.eval = @(s,t) chnk.helm2d.kern(zks(1), s, t, 'trep_diff', coefs(:,1)) - ...
                          chnk.helm2d.kern(zks(2), s, t, 'trep_diff', coefs(:,2));
        obj.fmm  = []; 
        if ( abs(coefs(2,1)-coefs(2,1)) < eps )
            obj.sing = 'log';
        else
            obj.sing = 'pv'; 
        end
    
    case {'trans_rep_prime','trep_p', 'trans_rep_p'}
        if ( nargin < 3 )
            coefs = ones(2,2);
        end
        obj.type = 'trep_p';
        obj.opdims = [1,2];
        obj.eval = @(s,t) chnk.helm2d.kern(zks(1), s, t, 'trep_p_diff', coefs(:,1)) - ...
                          chnk.helm2d.kern(zks(2), s, t, 'trep_p_diff', coefs(:,2));
        obj.fmm  = []; 
        if ( abs(coefs(2,1)-coefs(2,2)) < eps )
            obj.sing = 'log';
        else
            obj.sing = 'hs'; 
        end
    
    case {'c2trans', 'c2t', 'c2tr'}
        if ( nargin < 3 )
            coefs =  ones(2,2,2);
        end
        obj.type = 'c2tr';
        obj.opdims = [2,1];
        obj.eval = @(s,t) chnk.helm2d.kern(zks(1), s, t, 'c2tr_diff', coefs(:,:,1)) - ...
                          chnk.helm2d.kern(zks(2), s, t, 'c2tr_diff', coefs(:,:,2));
        obj.fmm  = []; 
        if ( abs(coefs(2,1,1)-coefs(2,1,2)) < eps )
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
