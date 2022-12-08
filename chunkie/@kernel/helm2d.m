function obj = helm2d(type, zk, eta)
%KERNEL.HELM2D   Construct the Helmholtz kernel.
%   KERNEL.HELM2D('s', ZK) or KERNEL.HELM2D('single', ZK) constructs the
%   single-layer Helmholtz kernel with wavenumber ZK.
%
%   KERNEL.HELM2D('d', ZK) or KERNEL.HELM2D('double', ZK) constructs the
%   double-layer Helmholtz kernel with wavenumber ZK.
%
%   KERNEL.HELM2D('sp', ZK) or KERNEL.HELM2D('sprime', ZK) constructs the
%   derivative of the single-layer Helmholtz kernel with wavenumber ZK.
%
%   KERNEL.HELM2D('c', ZK, ETA) or KERNEL.HELM2D('combined', ZK, ETA)
%   constructs the combined-layer Helmholtz kernel with wavenumber ZK and
%   parameter ETA, i.e., KERNEL.HELM2D('d', ZK) + 1i*ETA*KERNEL.HELM2D('s',
%   ZK).
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
obj.params.zk = zk;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 's');
        obj.fmm  = @(eps,s,t,sigma,pgt) chnk.helm2d.fmm(eps, zk, s, t, 's', sigma, pgt);
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'd');
        obj.fmm  = @(eps,s,t,sigma,pgt) chnk.helm2d.fmm(eps, zk, s, t, 'd', sigma, pgt);
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'sprime');
        % TODO: Add FMM for sprime.
        obj.sing = 'log';

    case {'c', 'combined'}
        if ( nargin < 3 )
            warning('Missing combined layer parameter eta. Defaulting to 1.');
            eta = 1;
        end
        obj.type = 'c';
        obj.params.eta = eta;
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'c', eta);
        obj.fmm  = @(eps,s,t,sigma,pgt) chnk.helm2d.fmm(eps, zk, s, t, 'c', sigma, pgt, eta);
        obj.sing = 'log';

    otherwise
        error('Unknown Helmholtz kernel type ''%s''.', type);

end

end
