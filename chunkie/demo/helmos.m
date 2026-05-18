function out = helmos(geo, bctype, zk, opts)
%HELMOS Helmholtz open-arc solver via Bruno-Lintner integral representation.
%
% Geometry-agnostic driver. Caller provides the chunkgraph (with all
% endpoints/corners/junctions as vertices) and a target list. Uses chunkie's
% RCIP machinery via chunkermat/chunkerkerneval with the open-arc identity
% override (opts.open_arc_eye=true) and rcipsav-aware field evaluation
% (opts.rcipsav).
%
% Status:
%   Dirichlet ('d')  -- delivers 10+ digits at typical grid distances
%                       (validated by helmos_smoke / helmos_subpanel_test).
%   Neumann   ('n')  -- formulation in place but limited to ~7 digits at
%                       present because chunkie's hypersingular 'dp'
%                       kernel does not yet have a kernel-split smooth
%                       quadrature like 'd'/'s' do; |sigma| grows with
%                       mesh refinement. Needs follow-up (port testhelmos
%                       Toper correction or add splitinfo to 'dp').
%
% Inputs:
%   geo - struct with fields:
%       .cg       - chunkgraph
%       .targets  - 2 x nt target locations
%       .ubdry_fn - function handle ubdry_fn(srcinfo, bctype) -> column,
%                   the boundary data for the BIE on the chunkgraph nodes
%   bctype - 'd' (Dirichlet) or 'n' (Neumann)
%   zk - wavenumber
%   opts - optional struct
%       .nsub        - RCIP refinement (default 30)
%       .tol         - GMRES tol (default 1e-12)
%       .gmres_maxit - default 200
%       .verbose     - print progress (default true)
%
% Output: struct with fields .sol, .iter, .ufield, .timings, .rcipsav,
%         .Keval. The Keval is returned so the caller can postprocess.

if nargin < 4 || isempty(opts), opts = []; end
nsub  = getfielddef(opts,'nsub',30);
tol   = getfielddef(opts,'tol',1e-12);
maxit = getfielddef(opts,'gmres_maxit',200);
vb    = getfielddef(opts,'verbose',true);
rcip_adap = getfielddef(opts,'rcip_adaptive_correction',false);

cg = geo.cg;
S  = kernel('helm','s',  zk);
T  = kernel('helm','dp', zk);
D  = kernel('helm','d',  zk);
Z  = kernel.zeros();
c  = 2.0;

% Bruno-Lintner formulation:
%   Dirichlet:   u = c S (-c T) sigma  -- field eval kernel = c S applied to tau
%   Neumann:     u = c D (-c S) sigma  -- field eval kernel = c D applied to tau
% Block matrix K + open_arc_eye gives [0 cA; cB I] form, where (A,B) =
% (S, T) for Dirichlet and (T, S) for Neumann.  After GMRES, the second
% density component is the auxiliary tau used by the field eval.
if bctype == 'd'
    K            = kernel([Z, c*S; c*T, Z]);
    Keval        = kernel([Z, c*S]);
    eval_type    = 's';        % chnk.kernsplit type for close eval
    fwlchs_layout = {{'zero',0}, {'s', c}; {'dp', c}, {'zero',0}};
elseif bctype == 'n'
    K            = kernel([Z, c*T; c*S, Z]);
    Keval        = kernel([Z, c*D]);   % NOTE: D, not T -- field uses double layer
    eval_type    = 'd';
    fwlchs_layout = {{'zero',0}, {'dp', c}; {'s', c}, {'zero',0}};
else
    error('HELMOS: bctype must be ''d'' or ''n''');
end

src = struct('r', cg.r(:,:), 'n', cg.n(:,:));
ubdry = geo.ubdry_fn(src, bctype, zk);

npts  = cg.npt;
nsys  = K.opdims(1)*npts;
rhs   = zeros(nsys,1);
rhs(1:K.opdims(1):end) = ubdry;

% forcewlchs covers all singular self+adj entries for the helmos block
% kernel, so GGQ self/adj corrections in chunkermat would be discarded.
% Use 'native' (smooth GL) and let wlchs+open_arc_eye fill in the
% diagonal/near-singular pieces.
matopts = struct('rcip',true,'nsub_or_tol',nsub, ...
                 'open_arc_eye',true,'rcip_savedepth',inf, ...
                 'rcip_adaptive_correction',rcip_adap, ...
                 'quad','native', ...
                 'forcewlchs', struct('zk',zk,'layout',{fwlchs_layout}));
matfree = getfielddef(opts,'matrix_free',false);
if matfree
    matopts.corrections = true;   % return only sparse cormat
end
if vb, fprintf('[helmos] building %s system (npts=%d, nsys=%d)...\n', ...
        ternary(matfree,'sparse','dense'), npts, nsys); end
t = tic; [A_or_cor, ~, rcipsav] = chunkermat(cg, K, matopts);
t_mat = toc(t);
if vb, fprintf('[helmos] matrix build: %.2f s\n', t_mat); end

if matfree
    chnkrmerge = merge(cg.echnks);
    src.r  = cg.r(:,:);  src.d  = cg.d(:,:);
    src.d2 = cg.d2(:,:); src.n  = cg.n(:,:);
    opdims_mat = [K.opdims(1); K.opdims(2)];
    fmm_smooth_opts = struct('forcefmm', getfielddef(opts,'forcefmm',true));
    apply = @(x) x + A_or_cor*x + ...
        chnk.chunkerkerneval_smooth(chnkrmerge, K, opdims_mat, x, src, [], fmm_smooth_opts);
else
    A = A_or_cor + eye(nsys);
    apply = @(x) A*x;
end

if vb, fprintf('[helmos] gmres (%s)...\n', ternary(matfree,'matrix-free','dense')); end
t = tic; [sol, flag, relres, iter] = gmres(apply, rhs, [], tol, maxit);
t_solve = toc(t);
if vb
    if iscell(iter), it = iter{end}; else, it = iter(end); end
    fprintf('[helmos] gmres: flag=%d relres=%.2e iter=%d in %.2f s\n', flag, relres, it, t_solve);
end

if vb, fprintf('[helmos] field evaluation at %d targets...\n', size(geo.targets,2)); end
t = tic;
% Single chunkerkerneval call: rcipsav + forcewlchs (block-kernel struct
% form) handles both the smooth far-field and the per-vertex per-edge
% kernel-split corner correction.  The layout selects which density
% component is the active one (column 2 = tau) and which kernsplit
% kernel type to apply (eval_type, with prefactor c).
eval_layout = {{'zero',0}, {eval_type, c}};
eval_opts = struct('rcipsav',  {rcipsav}, ...
                   'forcewlchs', struct('zk',zk,'layout',{eval_layout}), ...
                   'forcefmm', getfielddef(opts,'forcefmm',true));
ufield = chunkerkerneval(cg, Keval, sol, geo.targets, eval_opts);
t_eval = toc(t);
if vb, fprintf('[helmos] eval: %.2f s\n', t_eval); end

out = struct('sol',sol, 'iter',iter, 'ufield',ufield, ...
             'rcipsav',{rcipsav}, 'Keval',Keval, ...
             'timings',struct('matrix',t_mat,'solve',t_solve,'eval',t_eval));
end

function v = getfielddef(s,f,d)
    if isfield(s,f), v = s.(f); else, v = d; end
end

function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end
