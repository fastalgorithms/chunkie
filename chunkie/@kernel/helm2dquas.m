function obj = helm2dquas(type, zk, kappa, d, coefs, quad_opts)
%KERNEL.HELM2DQUAS   Construct the quasi periodic Helmholtz kernel.
%   KERNEL.HELM2DQUAS('s', ZK, KAPPA, D) or 
%   KERNEL.HELM2DQUAS('single', ZK, KAPPA, D) constructs the
%   single-layer quasiperiodic Helmholtz kernel with wavenumber ZK with 
%   period D and quasiperiodic factor kappa.
%
%   KERNEL.HELM2DQUAS('d', ZK, KAPPA, D) or 
%   KERNEL.HELM2DQUAS('double', ZK, KAPPA, D) constructs the
%   double-layer quasiperiodic Helmholtz kernel.
%
%   KERNEL.HELM2DQUAS('sp', ZK, KAPPA, D) or 
%   KERNEL.HELM2DQUAS('sprime', ZK, KAPPA, D) constructs the
%   derivative of the single-layer Helmholtz kernel.
%
%   KERNEL.HELM2DQUAS('c', ZK, KAPPA, D, COEFS) or
%   KERNEL.HELM2DQUAS('combined', ZK, KAPPA, D, COEFS)
%   constructs the combined-layer Helmholtz kernel with
%   parameter ETA, i.e., COEFS(1)*KERNEL.HELM2DQUAS('d', ZK) + 
%       COEFS(2)*KERNEL.HELM2DQUAS('s', ZK).
%
%   KERNEL.HELM2DQUAS('cp', ZK, KAPPA, D, COEFS) or 
%   KERNEL.HELM2DQUAS('combined', ZK, KAPPA, D, COEFS)
%   constructs the derivative of the combined-layer Helmholtz kernel with 
%   parameter ETA, i.e., COEFS(1)*KERNEL.HELM2DQUAS('dp', ZK, KAPPA, D) + 
%      COEFS(2)*KERNEL.HELM2DQUAS('sp', ZK, KAPPA, D).
%
%   if kappa is an array of values of length nkappa then the kernel 
%   becomes an (nkappa x 1) vector-valued kernel containing the scalar values
%   for each kappa. this is much more efficient than separate kernel 
%   evaluations for each kappa.
%   
%   all versions accept quad_opts to change the parameters in the lattice
%   sum computations. (see chnk.helm2dquas.latticecoefs)
%       quad_opts.l - (2) radius of periodic copies to exclude from 
%           lattice sum (excludes copies within l*d)
%       quad_opts.N - (40) number of lattice sums
%       quad_opts.m - (1e4) number of trapezoid rule quadrature nodes
%       quad_opts.a - (a) length of trpezoid quadrature rule
%   
% See also CHNK.HELM2DQUAS.KERN.

% author: Tristan Goodwill
  
if ( nargin < 1 )
    error('Missing Helmholtz kernel type.');
end

if ( nargin < 2 )
    error('Missing Helmholtz wavenumber.');
end

obj = kernel();
obj.name = 'quasiperiodic helmholtz';
obj.params.zk = zk;
obj.opdims = [numel(kappa) 1];


% lattice sum parameters
l=2; N = 40; a = 15; M = 1e4;

if nargin == 6
    if isfield(quad_opts,'l')
        l = quad_opts.l;
    end
    if isfield(quad_opts,'N')
        N = quad_opts.N;
    end
    if isfield(quad_opts,'M')
        M = quad_opts.M;
    end
    if isfield(quad_opts,'a')
        a = quad_opts.a;
    end
end

ns = (0:N).';
sn = chnk.helm2dquas.latticecoefs(ns,zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);

quas_param = [];
quas_param.kappa = kappa;
quas_param.d = d;
quas_param.l = l;
quas_param.sn = sn;

obj.params.quas_param = quas_param;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.helm2dquas.kern(zk, s, t, 's',quas_param);
        obj.fmm = [];
        obj.sing = 'log';
        if isscalar(kappa)
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0]};
        obj.splitinfo.action = {'r','r'};
        obj.splitinfo.functions = @(s,t) helm2dquas_s_split(zk,s,t,quas_param);
        end

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.helm2dquas.kern(zk, s, t, 'd',quas_param);
        obj.fmm = [];
        obj.sing = 'log';
        if isscalar(kappa)
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','r'};
        obj.splitinfo.functions = @(s,t) helm2dquas_d_split(zk,s,t,quas_param);
        end

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.helm2dquas.kern(zk, s, t, 'sprime',quas_param);
        obj.fmm = [];
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.helm2dquas.kern(zk, s, t, 'dprime',quas_param);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'c', 'combined'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'c';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2dquas.kern(zk, s, t, 'c',quas_param, coefs);
        obj.fmm = [];
        obj.sing = 'log';
        if isscalar(kappa)
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','r'};
        obj.splitinfo.functions = @(s,t) helm2dquas_c_split(zk,s,t,quas_param,coefs);
        end

    case {'cp', 'cprime'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'cprime';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2dquas.kern(zk, s, t, 'cprime',quas_param, coefs);
        obj.fmm = [];
        obj.sing = 'hs';

    otherwise
        error('Unknown Helmholtz kernel type ''%s''.', type);

end
end

function f = helm2dquas_s_split(zk,s,t,quas_param)
seval = chnk.helm2dquas.kern(zk, s, t, 's',quas_param);
s0eval = chnk.helm2d.kern(zk, s, t, 's');
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
logeval = log(abs(dist));
f = cell(2,1);
f{1} = seval + 1/(2*pi)*4*logeval.*imag(s0eval);
f{2} = 4*imag(s0eval);
end

function f = helm2dquas_d_split(zk,s,t,quas_param)
deval = chnk.helm2dquas.kern(zk, s, t, 'd',quas_param);
d0eval = chnk.helm2d.kern(zk, s, t, 'd');
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
logeval = log(abs(dist));
cauchyeval = (s.n(1,:)+1i*s.n(2,:))./(dist);
f = cell(3, 1);
f{1} = deval + 1/(2*pi)*4*logeval.*imag(d0eval) + 1/(2*pi)*real(cauchyeval);
f{2} = 4*imag(d0eval);
f{3} = ones(size(t.r,2),size(s.r,2));
end

function f = helm2dquas_c_split(zk,s,t,quas_param,coefs)
seval = chnk.helm2dquas.kern(zk, s, t, 's',quas_param);
deval = chnk.helm2dquas.kern(zk, s, t, 'd',quas_param);

s0eval = chnk.helm2d.kern(zk, s, t, 's');
d0eval = chnk.helm2d.kern(zk, s, t, 'd');
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
logeval = log(abs(dist));
cauchyeval = (s.n(1,:)+1i*s.n(2,:))./(dist);
f = cell(3, 1);
f{1} = coefs(1) * (deval + 2/pi*logeval.*imag(d0eval) + 1/(2*pi)*real(cauchyeval)) + ...
       coefs(2) * (seval + 2/pi*logeval.*imag(s0eval));
f{2} = coefs(1) * (4*imag(d0eval)) + coefs(2) * (4*imag(s0eval));
f{3} = coefs(1)*ones(size(t.r,2),size(s.r,2));
end

