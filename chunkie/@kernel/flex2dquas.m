function obj = flex2dquas(type, zk, kappa, d, coefs, quad_opts, ising)
%KERNEL.FLEX2DQUAS   Construct the quasi periodic flexural wave kernel.
%   KERNEL.FLEX2DQUAS('s', ZK, KAPPA, D) or
%   KERNEL.FLEX2DQUAS('single', ZK, KAPPA, D) constructs the
%   single-layer quasiperiodic flexural wave kernel with wavenumber ZK,
%   period D, and quasiperiodic factor KAPPA.
%
%   KERNEL.FLEX2DQUAS('d', ZK, KAPPA, D) or
%   KERNEL.FLEX2DQUAS('double', ZK, KAPPA, D) constructs the
%   double-layer quasiperiodic flexural wave kernel.
%
%   KERNEL.FLEX2DQUAS('sp', ZK, KAPPA, D) or
%   KERNEL.FLEX2DQUAS('sprime', ZK, KAPPA, D) constructs the
%   normal derivative of the single-layer quasiperiodic flexural wave kernel.
%
%   KERNEL.FLEX2DQUAS('clamped_plate_bcs', ZK, KAPPA, D) constructs the
%   boundary conditions applied to a point source for the clamped plate.
%
%   KERNEL.FLEX2DQUAS('clamped_plate_bcs_trx', ZK, KAPPA, D) constructs the
%   transmission boundary conditions to object for the clamped plate.
%
%   KERNEL.FLEX2DQUAS('clamped_plate', ZK, KAPPA, D) constructs the
%   kernels for the clamped plate integral equation.
%
%   KERNEL.FLEX2DQUAS('clamped_plate_eval', ZK, KAPPA, D) constructs the
%   clamped plate kernels for plotting.
%
%   KERNEL.FLEX2DQUAS('clamped_plate_eval_trx', ZK, KAPPA, D) constructs the
%   clamped plate kernels for plotting (transmission representation).
%
%   KERNEL.FLEX2DQUAS('free_plate_bcs', ZK, KAPPA, D, COEFS) constructs the
%   boundary conditions applied to a point source for the free plate.
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('free_plate_bcs_trx', ZK, KAPPA, D, COEFS) constructs
%   the transmission boundary conditions to object for the free plate.
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('free_plate', ZK, KAPPA, D, COEFS) constructs the
%   kernels for the free plate integral equation.
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('free_plate_eval', ZK, KAPPA, D, COEFS) constructs the
%   free plate kernels for plotting.
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('free_plate_eval_trx', ZK, KAPPA, D, COEFS) constructs
%   the free plate kernels for plotting (transmission representation).
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('supported_plate_bcs', ZK, KAPPA, D, COEFS) constructs
%   the boundary conditions applied to a point source for the supported plate.
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('supported_plate_bcs_trx', ZK, KAPPA, D, COEFS)
%   constructs the transmission boundary conditions to object for the
%   supported plate. COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('supported_plate', ZK, KAPPA, D, COEFS) constructs the
%   kernels for the supported plate integral equation.
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('supported_plate_eval', ZK, KAPPA, D, COEFS) constructs
%   the supported plate kernels for plotting.
%   COEFS(1) is the Poisson ratio nu.
%
%   KERNEL.FLEX2DQUAS('supported_plate_eval_trx', ZK, KAPPA, D, COEFS)
%   constructs the supported plate kernels for plotting (transmission
%   representation). COEFS(1) is the Poisson ratio nu.
%
%   if kappa is an array of values of length nkappa then the kernel
%   becomes an (nkappa x 1) vector-valued kernel containing the scalar
%   values for each kappa. this is much more efficient than separate
%   kernel evaluations for each kappa.
%
%   all versions accept quad_opts to change the parameters in the lattice
%   sum computations. (see chnk.helm2dquas.latticecoefs)
%       quad_opts.l - (2) radius of periodic copies to exclude from
%           lattice sum (excludes copies within l*d)
%       quad_opts.N - (40) number of lattice sums
%       quad_opts.m - (1e4) number of trapezoid rule quadrature nodes
%       quad_opts.a - (15) length of trapezoid quadrature rule
%
% See also CHNK.FLEX2DQUAS.KERN.

% author: Tristan Goodwill

if ( nargin < 1 )
    error('Missing flexural wave kernel type.');
end

if ( nargin < 2 )
    error('Missing flexural wave wavenumber.');
end

obj = kernel();
obj.name = 'quasiperiodic flexural wave';
obj.params.zk = zk;
obj.opdims = [numel(kappa) 1];

% lattice sum parameters
l = 2; N = 40; a = 15; M = 1e4;

if nargin >= 6
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

if nargin < 7
    ising = 1;
end

ns = (0:N).';
alpha = exp(1i*kappa*d);
Sn = chnk.flex2dquas.latticecoefs(ns,zk,d,kappa,alpha,a,M,l+1);

% lap2dquas lattice coefs (needed for free_plate kernels)
ns_l = (1:3).';
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs(ns_l,d,kappa,l+1);

obj.params.kappa = kappa;
obj.params.d = d;
obj.params.l = l;
obj.params.Sn = Sn;
obj.params.s0_l = s0_l;
obj.params.sn_l = sn_l;
obj.params.ising = ising;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 's', kappa, d, Sn, s0_l, sn_l, l, ising);
        obj.fmm = [];
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'sp', kappa, d, Sn, s0_l, sn_l, l, ising);
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'd', kappa, d, Sn, s0_l, sn_l, l, ising);
        obj.fmm = [];
        obj.sing = 'log';

    case {'clamped_plate_bcs'}
        obj.type = 'clamped_plate_bcs';
        obj.opdims = [2*numel(kappa) 1];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_bcs', kappa, d, Sn, s0_l, sn_l, l, ising);
        obj.fmm = [];
        obj.sing = 'log';

    case {'clamped_plate_bcs_trx'}
        obj.type = 'clamped_plate_bcs_trx';
        obj.opdims = [2*numel(kappa) 4];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_bcs_trx', kappa, d, Sn, s0_l, sn_l, l, ising);
        obj.fmm = [];
        obj.sing = 'log';

    case {'clamped_plate'}
        if ( nargin < 5 )
            coefs = [];
        end
        obj.type = 'clamped_plate';
        obj.params.coefs = coefs;
        obj.opdims = [2*numel(kappa) 2];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate', kappa, d, Sn, s0_l, sn_l, l, ising, coefs);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'clamped_plate_eval'}
        obj.type = 'clamped_plate_eval';
        obj.opdims = [numel(kappa) 2];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval', kappa, d, Sn, s0_l, sn_l, l, ising);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'clamped_plate_eval_trx'}
        obj.type = 'clamped_plate_eval_trx';
        obj.opdims = [4*numel(kappa) 2];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval_trx', kappa, d, Sn, s0_l, sn_l, l, ising);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'free_plate_bcs'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for free_plate_bcs.');
        end
        obj.type = 'free_plate_bcs';
        obj.params.coefs = coefs;
        obj.opdims = [2*numel(kappa) 1];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs', kappa, d, Sn, s0_l, sn_l, l, ising, coefs(1));
        obj.fmm = [];
        obj.sing = 'log';

    case {'free_plate_bcs_trx'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for free_plate_bcs_trx.');
        end
        obj.type = 'free_plate_bcs_trx';
        obj.params.coefs = coefs;
        obj.opdims = [2*numel(kappa) 4];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs_trx', kappa, d, Sn, s0_l, sn_l, l, ising, coefs(1));
        obj.fmm = [];
        obj.sing = 'log';

    case {'free_plate'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for free_plate.');
        end
        obj.type = 'free_plate';
        obj.params.coefs = coefs;
        obj.opdims = [4*numel(kappa) 2];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate', kappa, d, Sn, s0_l, sn_l, l, ising, coefs(1));
        obj.fmm = [];
        obj.sing = 'hs';

    case {'free_plate_eval'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for free_plate_eval.');
        end
        obj.type = 'free_plate_eval';
        obj.params.coefs = coefs;
        obj.opdims = [numel(kappa) 3];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval', kappa, d, Sn, s0_l, sn_l, l, ising, coefs(1));
        obj.fmm = [];
        obj.sing = 'hs';

    case {'free_plate_eval_trx'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for free_plate_eval_trx.');
        end
        obj.type = 'free_plate_eval_trx';
        obj.params.coefs = coefs;
        obj.opdims = [4*numel(kappa) 3];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval_trx', kappa, d, Sn, s0_l, sn_l, l, ising, coefs(1));
        obj.fmm = [];
        obj.sing = 'hs';

    case {'supported_plate_bcs'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for supported_plate_bcs.');
        end
        obj.type = 'supported_plate_bcs';
        obj.params.coefs = coefs;
        obj.opdims = [2*numel(kappa) 1];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_bcs', kappa, d, Sn, s0_l, sn_l, l, ising, coefs(1));
        obj.fmm = [];
        obj.sing = 'log';

    case {'supported_plate_bcs_trx'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for supported_plate_bcs_trx.');
        end
        obj.type = 'supported_plate_bcs_trx';
        obj.params.coefs = coefs;
        obj.opdims = [2*numel(kappa) 4];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_bcs_trx', kappa, d, Sn, s0_l, sn_l, l, ising, coefs(1));
        obj.fmm = [];
        obj.sing = 'log';

    case {'supported_plate'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for supported_plate.');
        end
        obj.type = 'supported_plate';
        obj.params.coefs = coefs;
        obj.opdims = [2*numel(kappa) 2];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate', kappa, d, Sn, s0_l, sn_l, l, ising, coefs);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'supported_plate_eval'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for supported_plate_eval.');
        end
        obj.type = 'supported_plate_eval';
        obj.params.coefs = coefs;
        obj.opdims = [numel(kappa) 2];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval', kappa, d, Sn, s0_l, sn_l, l, ising, coefs);
        obj.fmm = [];
        obj.sing = 'hs';

    case {'supported_plate_eval_trx'}
        if ( nargin < 5 )
            error('Missing Poisson ratio coefs for supported_plate_eval_trx.');
        end
        obj.type = 'supported_plate_eval_trx';
        obj.params.coefs = coefs;
        obj.opdims = [4*numel(kappa) 2];
        obj.eval = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval_trx', kappa, d, Sn, s0_l, sn_l, l, ising, coefs);
        obj.fmm = [];
        obj.sing = 'hs';

    otherwise
        error('Unknown quasiperiodic flexural wave kernel type ''%s''.', type);

end

end
