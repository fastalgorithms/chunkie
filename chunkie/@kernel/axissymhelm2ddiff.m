function obj = axissymhelm2ddiff(type, zks, coefs)
%KERNEL.AXISSYMHELM2DDIFF   Construct the difference of 
% axissymmetric Helmholtz kernels.
%   KERNEL.AXISSYMHELM2DDIFF('dp', ZKS) or 
%   KERNEL.AXISSYMHELM2DDIFF('dprime', ZKS) 
%   constructs the difference of Neumann data for double-layer 
%   axissymmetric Helmholtz kernels with wavenumbers ZKS(1) and ZKS(2), 
%   scaled by COEFS, i.e.
%   COEFS(1)*KERNEL.AXISSYMHELM2D('dp', ZKS(1)) -
%   COEFS(2)*KERNEL.AXISSYMHELM2D('dp', ZKS(2))
%
%  COEFS is optional. If not given then COEFS is set to [1 1].
%
%  NOTES: currently only supports coefs = [c c], and zks(1), real, and
%   zks(2) = I*zks(1)
%
% See also CHNK.AXISSYMHELM2D.KERN.

if ( nargin < 1 )
    error('Missing Helmholtz kernel type.');
end

if ( nargin < 2 )
    error('Missing Helmholtz wavenumber.');
end

zr1 = real(zks(1)); zi1 = imag(zks(1));
zr2 = real(zks(2)); zi2 = imag(zks(2));

ifreal = 0;
if zi1 ~=0 || zr2 ~=0 || zr1 ~= zi2
    ifreal = 1;
    if (zi1~= 0 || zi2 ~= 0)
        error('Wave numbers must be of the form zk, 1i*zk with zk real');
    end
end

obj = kernel();
obj.name = 'axissymhelmholtzdiff';
obj.params.zks = zks;
obj.opdims = [1 1];
if(nargin < 3)
    coefs = [1 1];
end


obj.params.coefs = coefs;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) coefs(1)*chnk.axissymhelm2d.kern(zks(1), ... 
                                    s, t, [0,0], 's') - ...
                          coefs(2)*chnk.axissymhelm2d.kern(zks(2), ... 
                                    s, t, [0,0], 's');
        obj.shifted_eval = @(s,t,o) coefs(1)*chnk.axissymhelm2d.kern(zks(1), ... 
                                    s, t, o, 's') - ...
                          coefs(2)*chnk.axissymhelm2d.kern(zks(2), ... 
                                    s, t, o, 's');
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) coefs(1)*chnk.axissymhelm2d.kern(zks(1), ... 
                                    s, t, [0,0], 'd') - ...
                          coefs(2)*chnk.axissymhelm2d.kern(zks(2), ... 
                                    s, t, [0,0], 'd');
        obj.shifted_eval = @(s,t,o) coefs(1)*chnk.axissymhelm2d.kern(zks(1), ... 
                                    s, t, o, 'd') - ...
                          coefs(2)*chnk.axissymhelm2d.kern(zks(2), ... 
                                    s, t, o, 'd');
        obj.fmm = [];
        obj.sing = 'log';       

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) coefs(1)*chnk.axissymhelm2d.kern(zks(1), ... 
                                    s, t, [0,0], 'sprime') - ...
                          coefs(2)*chnk.axissymhelm2d.kern(zks(2), ... 
                                    s, t, [0,0], 'sprime');
        obj.shifted_eval = @(s,t,o) coefs(1)*chnk.axissymhelm2d.kern(zks(1), ... 
                                    s, t, o, 'sprime') - ...
                          coefs(2)*chnk.axissymhelm2d.kern(zks(2), ... 
                                    s, t, o, 'sprime');
        obj.fmm = [];
        obj.sing = 'log';

    case {'dp', 'dprime'}

        if coefs(1) ~= coefs(2)
            error('Coefs must be of form [c c]');   
        end
        obj.type = 'dp';
        if(~ifreal)
        obj.eval = @(s,t) coefs(1)*chnk.axissymhelm2d.kern(real(zks(1)),  ... 
                                                s, t, [0,0], 'dprimediff');
        obj.shifted_eval = @(s,t,o) coefs(1)*chnk.axissymhelm2d.kern(real(zks(1)),  ... 
                                                s, t, o, 'dprimediff');
        obj.fmm = [];
        else
        obj.eval = @(s,t) coefs(1)*chnk.axissymhelm2d.kern(zks,  ... 
                                                s, t, [0,0], 'dprime_re_diff');
        obj.shifted_eval = @(s,t,o) coefs(1)*chnk.axissymhelm2d.kern(zks,  ... 
                                                s, t, o, 'dprime_re_diff');
        obj.fmm = [];
        end
        if ( abs(coefs(1)-coefs(2)) < eps )
            obj.sing = 'log';
        else
            obj.sing = 'hs';
        end

    otherwise
        error('Unknown axissymmetric Helmholtz difference kernel type ''%s''.', type);

end

end
