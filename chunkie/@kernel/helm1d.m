function obj = helm1d(type, zk, coefs)
%KERNEL.HELM2D   Construct the Helmholtz kernel.
%   KERNEL.HELM2D('s', ZK) or KERNEL.HELM2D('single', ZK) constructs the
%   single-layer Helmholtz kernel with wavenumber ZK.
%
%
% See also CHNK.HELM2D.KERN.

if ( nargin < 1 )
    error('Missing Helmholtz kernel type.');
end

if ( nargin < 2 )
    error('Missing Helmholtz wavenumber.');
end

obj = kernel();
obj.name = 'helmholtz1d';
obj.params.zk = zk;
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.helm1d.kern(zk, s, t, 's');
        obj.fmm  = []; 
        obj.sing = 'removable';

    otherwise
        error('Unknown Helmholtz kernel type ''%s''.', type);

end

end
