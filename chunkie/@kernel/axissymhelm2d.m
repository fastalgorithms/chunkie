function obj = axissymhelm2d(type, zk, coefs)
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
    error('Missing Helmholtz kernel type.');
end

if ( nargin < 2 )
    error('Missing Helmholtz wavenumber.');
end


zr = real(zk); zi = imag(zk);
if abs(zr)*abs(zi) > eps
    error('Only purely real or purely imaginary wavenumbers supported');
end

obj = kernel();
obj.name = 'axissymhelmholtz';
obj.params.zk = zk;
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.axissymhelm2d.kern(zk, s, t, [0,0], 's');
        obj.shifted_eval = @(s,t,o) chnk.axissymhelm2d.kern(zk, s, t, o, 's');
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.axissymhelm2d.kern(zk, s, t, [0,0], 'd');
        obj.shifted_eval = @(s,t,o) chnk.axissymhelm2d.kern(zk, s, t, o, 'd');
        obj.fmm = [];
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.axissymhelm2d.kern(zk, s, t, [0,0], 'sprime');
        obj.shifted_eval = @(s,t,o) chnk.axissymhelm2d.kern(zk, s, t, o, 'sprime');
        obj.fmm = [];
        obj.sing = 'log';

    case {'c', 'combined'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'c';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.axissymhelm2d.kern(zk, s, t, [0,0], 'c', coefs);
        obj.shifted_eval = @(s,t,o) chnk.axissymhelm2d.kern(zk, s, t, o, 'c', coefs);
        obj.fmm = [];
        obj.sing = 'log';
     case {'neu_rpcomb'}
        obj.type = 'neu_rpcomb';
        if ( nargin < 3 )
            warning('Missing coefficient of 1i*D. Defaulting to 1.');
            coefs = 1;
        end
        obj.eval = @(s,t) chnk.axissymhelm2d.kern(zk, s, t, [0,0], 'neu_rpcomb', coefs);
        obj.shifted_eval = @(s,t,o) chnk.axissymhelm2d.kern(zk, s, t, o, 'neu_rpcomb', coefs);
        obj.fmm = [];
        obj.sing = 'log';
        obj.opdims = [3 3];
        obj.params.c1 = -1.0/(0.5 + 0.25*1i*coefs);
        obj.params.c2 = -1i*coefs/(0.5 + 0.25*1i*coefs);

    otherwise
        error('Unknown axissym Helmholtz kernel type ''%s''.', type);

end

end
