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
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 's');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 's', sigma);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0]};
        obj.splitinfo.action = {'r','r'};
        obj.splitinfo.function = cell(2,1);
        obj.splitinfo.function{1} = @(s,t) obj.eval(s,t) ...
                                         + 1/(2*pi)*4*log(abs((s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)'))).*imag(obj.eval(s,t));
        obj.splitinfo.function{2} = @(s,t) 4*imag(obj.eval(s,t));

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'd');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'd', sigma);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','r'};
        obj.splitinfo.function = cell(3,1);
        obj.splitinfo.function{1} = @(s,t) obj.eval(s,t) ...
                                           + 1/(2*pi)*4*log(abs((s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)'))).*imag(obj.eval(s,t)) ... 
                                           + 1/(2*pi)*real((s.n(1,:)+1i*s.n(2,:))./((s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)')));
        obj.splitinfo.function{2} = @(s,t) 4*imag(obj.eval(s,t));
        obj.splitinfo.function{3} = @(s,t) ones(size(t.r,2),size(s.r,2));

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'sprime');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'sprime', sigma);
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'dprime');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'dprime', sigma);
        obj.sing = 'hs';

    case {'c', 'combined'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'c';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'c', coefs);
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'c', sigma, coefs);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','r'};
        obj.splitinfo.functions = @(s,t) helm2d_c(zk,s,t,coefs);

    case {'cp', 'cprime'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'cprime';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'cprime', coefs);
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'cprime', sigma, coefs);
        obj.sing = 'hs';

    otherwise
        error('Unknown Helmholtz kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end


end

function f = helm2d_c(zk,s,t,coefs)
seval = chnk.helm2d.kern(zk, s, t, 's');
deval = chnk.helm2d.kern(zk, s, t, 'd');
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
logeval = log(abs(dist));
cauchyeval = (s.n(1,:)+1i*s.n(2,:))./(dist);
f = cell(3, 1);
f{1} = coefs(1) * (deval + 2/pi*logeval.*imag(deval) + 1/(2*pi)*real(cauchyeval)) + ...
       coefs(2) * (seval + 2/pi*logeval.*imag(seval));
f{2} = coefs(1) * (4*imag(deval)) + coefs(2) * (4*imag(seval));
f{3} = coefs(1)*ones(size(t.r,2),size(s.r,2));
end

