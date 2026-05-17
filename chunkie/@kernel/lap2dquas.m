function obj = lap2dquas(type, kappa, d, coefs, quad_opts)
%KERNEL.LAP2DQUAS   Construct the quasi periodic Laplace kernel.
%   KERNEL.LAP2DQUAS('s', KAPPA, D) or
%   KERNEL.LAP2DQUAS('single', KAPPA, D) constructs the
%   single-layer quasiperiodic Laplace kernel with period D and
%   quasiperiodic factor KAPPA.
%
%   KERNEL.LAP2DQUAS('d', KAPPA, D) or
%   KERNEL.LAP2DQUAS('double', KAPPA, D) constructs the
%   double-layer quasiperiodic Laplace kernel.
%
%   KERNEL.LAP2DQUAS('sp', KAPPA, D) or
%   KERNEL.LAP2DQUAS('sprime', KAPPA, D) constructs the
%   normal derivative of the single-layer quasiperiodic Laplace kernel.
%
%   KERNEL.LAP2DQUAS('st', KAPPA, D) or
%   KERNEL.LAP2DQUAS('stau', KAPPA, D) constructs the
%   tangential derivative of the single-layer quasiperiodic Laplace kernel.
%
%   KERNEL.LAP2DQUAS('hilb', KAPPA, D) constructs the
%   quasiperiodic Hilbert transform kernel.
%
%   KERNEL.LAP2DQUAS('hilbprime', KAPPA, D) constructs the
%   tangential derivative of the quasiperiodic Hilbert transform kernel.
%
%   KERNEL.LAP2DQUAS('dp', KAPPA, D) or
%   KERNEL.LAP2DQUAS('dprime', KAPPA, D) constructs the
%   normal derivative of the double-layer quasiperiodic Laplace kernel.
%
%   if kappa is an array of values of length nkappa then the kernel
%   becomes an (nkappa x 1) vector-valued kernel containing the scalar
%   values for each kappa. this is much more efficient than separate
%   kernel evaluations for each kappa.
%
%   all versions accept quad_opts to change the parameters in the lattice
%   sum computations. (see chnk.lap2dquas.latticecoefs)
%       quad_opts.l - (2) radius of periodic copies to exclude from
%           lattice sum (excludes copies within l*d)
%
% See also CHNK.LAP2DQUAS.KERN.

% author: Tristan Goodwill

if ( nargin < 1 )
    error('Missing Laplace kernel type.');
end

if ( nargin < 2 )
    error('Missing quasiperiodic parameter kappa.');
end

obj = kernel();
obj.name = 'quasiperiodic laplace';
obj.opdims = [numel(kappa) 1];

% lattice sum parameters
l = 2; N = 40;
if nargin >= 5
    if isfield(quad_opts,'l')
        l = quad_opts.l;
    end
    if isfield(quad_opts,'N')
        N = quad_opts.N;
    end
end

ns = (1:N);
[s0,sn] = chnk.lap2dquas.latticecoefs(ns,d,kappa,l);

quas_param = [];
quas_param.kappa = kappa;
quas_param.d = d;
quas_param.l = l;
quas_param.s0 = s0;
quas_param.sn = sn;

obj.params.quas_param = quas_param;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.lap2dquas.kern(s, t, 's', kappa, d, s0, sn, l, 1);
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.lap2dquas.kern(s, t, 'd', kappa, d, s0, sn, l, 1);
        obj.fmm = [];
        obj.sing = 'smooth';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.lap2dquas.kern(s, t, 'sp', kappa, d, s0, sn, l, 1);
        obj.fmm = [];
        obj.sing = 'smooth';

    case {'st', 'stau'}
        obj.type = 'st';
        obj.eval = @(s,t) chnk.lap2dquas.kern(s, t, 'stau', kappa, d, s0, sn, l, 1);
        obj.fmm = [];
        obj.sing = 'pv';

    case {'hilb'}
        obj.type = 'hilb';
        obj.eval = @(s,t) chnk.lap2dquas.kern(s, t, 'hilb', kappa, d, s0, sn, l, 1);
        obj.fmm = [];
        obj.sing = 'pv';

    case {'hilbprime'}
        obj.type = 'hilbprime';
        obj.eval = @(s,t) chnk.lap2dquas.kern(s, t, 'hilbprime', kappa, d, s0, sn, l, 1);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.lap2dquas.kern(s, t, 'dp', kappa, d, s0, sn, l, 1);
        obj.fmm = [];
        obj.sing = 'hs';

    otherwise
        error('Unknown quasiperiodic Laplace kernel type ''%s''.', type);

end

end
