function obj = helm2d(type, zk, coefs)
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
%   KERNEL.HELM2D('c', ZK, COEFS) or KERNEL.HELM2D('combined', ZK, COEFS)
%   constructs the combined-layer Helmholtz kernel with wavenumber ZK and
%   parameter ETA, i.e., COEFS(1)*KERNEL.HELM2D('d', ZK) + 
%   COEFS(2)*KERNEL.HELM2D('s', ZK).
%
%   KERNEL.HELM2D('cp', ZK, COEFS) or KERNEL.HELM2D('combined', ZK, COEFS)
%   constructs the derivative of the combined-layer Helmholtz kernel with 
%   wavenumber ZK and  parameter ETA, i.e., COEFS(1)*KERNEL.HELM2D('dp', ZK) + 
%   COEFS(2)*KERNEL.HELM2D('sp', ZK).
% See also CHNK.HELM2D.KERN.

% author: Dan Fortunato
  
if ( nargin < 1 )
    error('Missing Helmholtz kernel type.');
end

if ( nargin < 2 )
    error('Missing Helmholtz wavenumber.');
end

obj = kernel();
obj.name = 'helmholtz';
obj.params.zk = zk;
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.helm2d_sheet.kern(zk, s, t, 's');
        obj.fmm  = [];
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0]};
        obj.splitinfo.action = {'r','r'};
        obj.splitinfo.functions = @(s,t) helm2d_s_split(zk,s,t);

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.helm2d_sheet.kern(zk, s, t, 'd');
        obj.fmm  = [];
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','r'};

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.helm2d_sheet.kern(zk, s, t, 'sprime');
        obj.fmm  = [];
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.helm2d_sheet.kern(zk, s, t, 'dprime');
        obj.sing = 'hs';

    case {'c', 'combined'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'c';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2d_sheet.kern(zk, s, t, 'c', coefs);
        obj.sing = 'log';
        obj.splitinfo = [];

    case {'cp', 'cprime'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'cprime';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2d_sheet.kern(zk, s, t, 'cprime', coefs);
        obj.sing = 'hs';

    otherwise
        error('Unknown Helmholtz kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end


end
