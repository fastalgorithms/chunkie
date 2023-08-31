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

if zi1 ~=0 || zr2 ~=0 || zr1 ~= zi2
    error('Wave numbers must be of the form zk, 1i*zk with zk real');
end

obj = kernel();
obj.name = 'axissymhelmholtzdiff';
obj.params.zks = zks;
obj.opdims = [1 1];
if(nargin < 3)
    coefs = [1 1];
end

if coefs(1) ~= coefs(2)
    error('Coefs must be of form [c c]');   
end

obj.params.coefs = coefs;

switch lower(type)


    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) coefs(1)*chnk.axissymhelm2d.kern(real(zks(1)),  ... 
                                                s, t, [0,0], 'dprimediff');
        obj.shifted_eval = @(s,t,o) coefs(1)*chnk.axissymhelm2d.kern(real(zks(1)),  ... 
                                                s, t, o, 'dprimediff');
        obj.fmm = [];
        if ( abs(coefs(1)-coefs(2)) < eps )
            obj.sing = 'log';
        else
            obj.sing = 'hs';
        end

    otherwise
        error('Unknown axissymmetric Helmholtz difference kernel type ''%s''.', type);

end

end
