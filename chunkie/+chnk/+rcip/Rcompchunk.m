function [R,rcipsav]=Rcompchunk(chnkr,iedgechunks,fkern,ndim,vert0, ...
    Pbc,PWbc,nsub,starL,circL,starS,circS,ilist,starL1,circL1,...
    sbclmat,sbcrmat,lvmat,rvmat,u,opts)
%CHNK.RCIP.Rcompchunk carry out the forward recursion for computing
% the preconditioner R where geometry is described as a chunker
%
% This routine is not intended to be user-callable 
%
% Function is passed as a handle, number of equations is given by
% ndim
%
% Kernel on input takes in arguments (chnkrlocal,ilistl);
% 
% Note that matrix must be scaled to have identity on the diagonal,
% will not work with scaled version of identity
%

% author: Shidong Jiang, derived from Rcomp.m by Johan Helsing in the
% RCIP tutorial 
% modified: Jeremy Hoskins, Manas Rachh


k = chnkr.k;  
dim = chnkr.dim;
nedge = size(iedgechunks,2);

glxs = chnkr.tstor;
glws = chnkr.wstor;

% return what's needed to interpolate from coarse

rcipsav = [];
rcipsav.k = k;
rcipsav.ndim = ndim;
rcipsav.nedge = nedge;
rcipsav.Pbc = Pbc;
rcipsav.PWbc = PWbc;
rcipsav.starL = starL;
rcipsav.starL1 = starL1;
rcipsav.starS = starS;
rcipsav.circL = circL;
rcipsav.circL1 = circL1;
rcipsav.circS = circS;
rcipsav.ilist = ilist;
rcipsav.nsub = nsub;

if (nargin < 16 || isempty(sbclmat) || isempty(sbcrmat) || ...
        isempty(lvmat) || isempty(rvmat) || isempty(u))
    [sbclmat,sbcrmat,lvmat,rvmat,u] = chnk.rcip.shiftedlegbasismats(k); 
end

if nargin < 21
    opts = [];
end

savedepth = 10;
if isfield(opts,'rcip_savedepth')
    savedepth = opts.rcip_savedepth;
end
savedepth = max(savedepth,0);
savedepth = min(savedepth,nsub);

% open_arc_eye Bruno-Lintner structure makes the (1,1) block of the local
% system matrix near-singular; use a Schur-complement-on-the-other-block
% inverse (chnk.rcip.myinv) inside SchurBana to avoid ~1e-23 RCOND warnings.
use_myinv = isfield(opts,'open_arc_eye') && opts.open_arc_eye;
% sub-block size for the myinv 2x2-block-of-(d x d) Schur structure
% (1 for helmos, 2 for open-arc Stokes mobility).
myinv_subdim = 1;
if isfield(opts,'open_arc_eye_subdim') && ~isempty(opts.open_arc_eye_subdim)
    myinv_subdim = opts.open_arc_eye_subdim;
end

rcipsav.savedepth = savedepth;
rcipsav.use_myinv = use_myinv;
rcipsav.myinv_subdim = myinv_subdim;

rcipsav.R = cell(nsub+1,1);
rcipsav.MAT = cell(nsub,1);
rcipsav.chnkrlocals = cell(nsub,1);

% Optional: cache full per-level MAT to enable srhs_recurse_only without
% rebuilding the (whole) chunkermat. Used by callers that re-run the
% singular-RHS recursion many times with different srhs_eval (e.g.,
% I2I solves with multiple data_in inputs sharing the same operator).
cache_mat_for_srhs = isfield(opts,'cache_mat_for_srhs') && opts.cache_mat_for_srhs;
if cache_mat_for_srhs
    rcipsav.MAT_full      = cell(nsub,1);
    rcipsav.chnkrlocals_array = cell(nsub,1);   % un-merged (array of nedge)
end

% grab only those kernels relevant to this vertex

if(size(fkern)==1)
    fkernlocal = fkern;
else
  fkernlocal(nedge,nedge) = kernel();
    for i=1:nedge
        ici = iedgechunks(1,i);
        for j=1:nedge
            icj = iedgechunks(1,j);
            fkernlocal(i,j) = fkern(ici,icj);
        end
    end

end

rcipsav.fkernlocal = fkernlocal;
rcipsav.iedgechunks = iedgechunks;   % needed for external srhs_recurse_only
rcipsav.vert0       = vert0;          % corner vertex (global coord)

% Singular-RHS RCIP support (Helsing & Karlsson, eq. 28). When opts.srhs_eval
% is provided, run the r_f^star vector recursion in tandem with the standard
% R recursion. The callback returns f at the 96*ndim GL nodes of the level-i
% type-b 6-panel local mesh (in global coordinates).
srhs_eval = [];
if isfield(opts,'srhs_eval') && ~isempty(opts.srhs_eval)
    srhs_eval = opts.srhs_eval;
end
rcipsav.has_srhs = ~isempty(srhs_eval);

% Debug instrumentation for the singular-RHS recursion. When set, prints
% per-level diagnostics and saves a per-level trace into rcipsav.srhs_trace.
srhs_debug = isfield(opts,'srhs_debug') && opts.srhs_debug;
if srhs_debug
    rcipsav.srhs_trace = {};
end

% get coefficients of recentered edge chunks and figure out orientation

km1 = k-1;
rcs = zeros(km1,dim,nedge);
dcs = zeros(k,dim,nedge);
d2cs = zeros(k,dim,nedge);
dscal = zeros(nedge,1);
d2scal = zeros(nedge,1);
ctr = zeros(dim,nedge);

ileftright = zeros(nedge,1);
nextchunk = zeros(nedge,1);

for i = 1:nedge
    ic = iedgechunks(1,i);
    ie = iedgechunks(2,i);
    chnkri = chnkr(ic);
    r = chnkri.r(:,:,ie);
    d = chnkri.d(:,:,ie);    
    d2 = chnkri.d2(:,:,ie);    
    il = chnkri.adj(1,ie);
    ir = chnkri.adj(2,ie);
    if (il > 0 && ir < 0)
        nextchunk(i) = il;
        ileftright(i) = 1;

        rr = rvmat*(r.'); r = r - rr(:);
        ctr(:,i) = rr;
        rcs(:,:,i) = sbcrmat*(r.');
        dcs(:,:,i) = u*(d.');
        d2cs(:,:,i) = u*(d2.');
        dscal(i) = 2;
        d2scal(i) = 4; 
    elseif (il < 0 && ir > 0)
        nextchunk(i) = ir;
        ileftright(i) = -1;
        rl = lvmat*(r.'); r = r - rl(:);
        ctr(:,i) = rl;
        rcs(:,:,i) = sbclmat*(r.');
        dcs(:,:,i) = u*(d.');
        d2cs(:,:,i) = u*(d2.');
        dscal(i) = 2;
        d2scal(i) = 4; 
    else
        error('RCIP: edge chunk not adjacent to one vertex and one neighbor')
    end
end

rcipsav.ctr = ctr;
rcipsav.rcs = rcs;
rcipsav.dcs = dcs;
rcipsav.d2cs = d2cs;
rcipsav.dscal = dscal;
rcipsav.d2scal = d2scal;
rcipsav.ileftright = ileftright;
rcipsav.glxs = glxs;
rcipsav.glws = glws;

pref = []; 
pref.k = k;
pref.nchstor = 5;

R = [];

% size of the system matrix
nsys = 3*k*nedge*ndim;
  
% size of the preconditioner R
nR = 2*k*nedge*ndim;

ts = cell(nedge,1);
pref = [];
pref.k = k;
chnkrlocal(1,nedge) = chunker(pref);



if(size(fkern)==1)
    fkernlocal = fkern;
    if isa(fkern, 'kernel')
        if isa(fkern.shifted_eval, 'function_handle')
            fkernlocal.eval = @(s,t) fkern.shifted_eval(s, t, vert0);
        end
    end
else
    fkernlocal(nedge,nedge) = kernel();
    for i=1:nedge
        ici = iedgechunks(1,i);
        for j=1:nedge
            icj = iedgechunks(1,j);
            fkernlocal(i,j) = fkern(ici,icj);
            if isa(fkern(ici,icj), 'kernel')
                if isa(fkern(ici,icj).shifted_eval, 'function_handle')
                    fkernlocal(i,j).eval = ...
                      @(s,t) fkern(ici,icj).shifted_eval(s,t,vert0);
                end
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin recursion proper

r_f_star = [];     % singular-RHS recursion state (only used when srhs_eval set)

h0=ones(nedge,1);
for level=1:nsub
    h = h0/2^(nsub-level);

    for i=1:nedge
        if ileftright(i) == -1
            if level == nsub 
                ts{i} =  [0, 0.5, 1]*h(i);
            else
                ts{i} =  [0, 0.5, 1, 2]*h(i);
            end
        else
            if level == nsub
                ts{i} = -[1, 0.5, 0]*h(i);
            else
                ts{i} = -[2, 1, 0.5, 0]*h(i);
            end
            
        end
    end
    % construct local chunks around the corner
    for i=1:nedge
        chnkrlocal(i) = chnk.rcip.chunkerfunclocal(@(t) shiftedcurve(t,rcs(:,:,i),dcs(:,:,i), ...
            dscal(i),d2cs(:,:,i),d2scal(i),ileftright(i)),ts{i},pref,glxs,glws);
    end
    
    % at the top level, append/prepend the next chunk
    if level == nsub
        for i = 1:nedge
            ic = iedgechunks(1,i);
            nc = nextchunk(i);
            nchi = chnkrlocal(i).nch;
            chnkrlocal(i) = chnkrlocal(i).addchunk(1);
            chnkrlocal(i).r(:,:,nchi+1) = chnkr(ic).r(:,:,nc)-ctr(:,i);
            chnkrlocal(i).d(:,:,nchi+1) = chnkr(ic).d(:,:,nc);                
            chnkrlocal(i).d2(:,:,nchi+1) = chnkr(ic).d2(:,:,nc); 

            if ileftright(i) == -1
                chnkrlocal(i).adj(1,nchi+1) = nchi;
                chnkrlocal(i).adj(2,nchi+1) = -1;
                chnkrlocal(i).adj(2,nchi) = nchi+1;
            else
                chnkrlocal(i).adj(1,nchi+1) = -1;
                chnkrlocal(i).adj(2,nchi+1) = 1;
                chnkrlocal(i).adj(1,1) = nchi+1;
                chnkrlocal(i) = chnkrlocal(i).sort();
            end
            chnkrlocal(i).n = normals(chnkrlocal(i));
            chnkrlocal(i).wts = weights(chnkrlocal(i));
        end
    end

    
% construct the system matrix for local chunks
    if level == 1
        ilistl = [];
    else
        ilistl = ilist;
    end

    % test for opdims ~= [1,1]
    [MAT,opts] = chunkermat(chnkrlocal, fkernlocal, opts, ilistl);
    % Apply cross-edge wLCHS at every recursive level (not just level==1).
    % Each refinement level has its own cross-edge near interactions across
    % the corner, all of which need the kernel-split correction.
    if nedge > 1
        MAT = add_level1_cross_wlchs(MAT, chnkrlocal, ndim, opts, fkernlocal);
    end
    

%
    MAT = eye(nsys) + MAT;
    if level==1    %  Dumb and lazy initializer for R, for now
  %R=eye(nR);
        if use_myinv
            R = chnk.rcip.myinv(MAT(starL,starL), myinv_subdim);
        else
            R = inv(MAT(starL,starL));
        end
        if level >= nsub-savedepth+1
            rcipsav.R{1} = R;
        end
    end
    if savedepth < nsub && level == nsub-savedepth+1
        rcipsav.R{level} = R;
    end

    % --- Singular-RHS RCIP forward recursion, vector form (Helsing 2022 eq 28) ---
    % Run before the SchurBana update of R, since Rfstep needs R_{i-1}.
    if ~isempty(srhs_eval)
        % Translate the local 6-panel mesh to global coords for the callback.
        chnkrlocal_global = chnkrlocal;
        for ie = 1:nedge
            chnkrlocal_global(ie) = chnkrlocal(ie) + ctr(:,ie);
        end
        b_ib = srhs_eval(chnkrlocal_global, vert0, iedgechunks(1,:), level, ndim);
        b_ib = b_ib(:);
        if numel(b_ib) ~= nsys
            error(['CHNK.RCIP.RCOMPCHUNK: srhs_eval returned %d entries, ' ...
                   'expected %d (3*k*nedge*ndim).'], numel(b_ib), nsys);
        end
        if level == 1
            % Lazy init: r_{f,0}^star = R_0 * b_{1b}(starL).
            r_f_star = R * b_ib(starL);
        end
        r_f_star_new = chnk.rcip.Rfstep(MAT, R, r_f_star, b_ib(circL), ...
            Pbc, PWbc, starL, circL, starS, circS, use_myinv);

        if srhs_debug
            tr.level = level;
            tr.norm_R       = norm(R);
            tr.norm_M       = norm(MAT);
            tr.norm_b_starL = norm(b_ib(starL));
            tr.norm_b_circL = norm(b_ib(circL));
            tr.norm_r_f_star = norm(r_f_star_new);
            tr.b_ib = b_ib;                     % full level-i b vector
            tr.r_f_star_new = r_f_star_new;
            tr.r = arrayfun(@(c) c.r, chnkrlocal_global, 'UniformOutput', false);
            rcipsav.srhs_trace{end+1} = tr; %#ok<AGROW>
        end
    end

    R=chnk.rcip.SchurBana(Pbc,PWbc,MAT,R,starL,circL,starS,circS,use_myinv,myinv_subdim);
    if level >= nsub-savedepth+1
        rcipsav.R{level+1} = R;
        rcipsav.MAT{level} = MAT(starL,circL);
        rcipsav.chnkrlocals{level} = merge(chnkrlocal);
        if cache_mat_for_srhs
            rcipsav.MAT_full{level}          = MAT;
            rcipsav.chnkrlocals_array{level} = chnkrlocal;
        end
    end

    if ~isempty(srhs_eval)
        r_f_star = r_f_star_new;
    end
end

if ~isempty(srhs_eval)
    rcipsav.r_f_star = r_f_star;
end

end

function [r,d,d2] = shiftedcurve(t,rc,dc,scald,d2c,scald2,ilr)
n = length(dc);
nm1 = n-1;
if (ilr == -1)
    tt = 2*t-1;
else
    tt = 2*t+1;
end
pols = lege.pols(tt,nm1); pols = pols.';
r = (t.*(pols(:,1:nm1)*rc)).';
d = (scald*(pols*dc)).';
d2 = (scald2*(pols*d2c)).';

end

function MAT = add_level1_cross_wlchs(MAT, chnkrlocal, ndim, opts, fkernlocal)
% Apply off-edge near (cross-edge) corrections at first RCIP level.
%
% Two modes supported:
%  - opts.forcewlchs is struct with explicit per-(r,c) layout: same layout
%    used for all edge pairs (legacy testhelmos / matrix-valued op path).
%  - opts.forcewlchs == true (boolean) AND fkernlocal is a kernel cell
%    array: per-(target_edge, source_edge) layout is auto-derived by
%    inspecting fkernlocal(it, is) (handles 'd', 's', 'sp', 'dp', 'c'/'combined').

if ~isfield(opts,'forcewlchs')
    return
end
fw = opts.forcewlchs;
have_struct = isstruct(fw) && isfield(fw,'layout') && isfield(fw,'zk');
have_bool   = (islogical(fw) || (isnumeric(fw) && isscalar(fw))) && fw && ...
              nargin >= 5 && ~isempty(fkernlocal);
if ~have_struct && ~have_bool
    return
end

nedge = numel(chnkrlocal);
npts = arrayfun(@(c) c.npt, chnkrlocal);
offset = [0, cumsum(npts(1:end-1))*ndim];

if have_struct
    layout = fw.layout;
    zk = fw.zk;
    if ~iscell(layout) || size(layout,1) ~= ndim || size(layout,2) ~= ndim
        return
    end
end

for it = 1:nedge
    for is = 1:nedge
        if it == is
            continue
        end
        for r = 1:ndim
            for c = 1:ndim
                if have_struct
                    entry = layout{r,c};
                    if isempty(entry); continue; end
                    if iscell(entry)
                        typ = entry{1}; coef = entry{2};
                    elseif isstruct(entry)
                        typ = entry.type; coef = entry.coef;
                    else
                        continue
                    end
                    if strcmpi(typ,'zero'); continue; end
                    [types, coefs] = decompose_kernel_spec(typ, coef);
                else
                    % Boolean mode: derive from fkernlocal(it, is). When
                    % fkernlocal is a single (scalar) kernel object, use it
                    % uniformly for every (it, is) pair.
                    if isscalar(fkernlocal)
                        bk = fkernlocal;
                    else
                        bk = fkernlocal(it, is);
                    end
                    if ~isa(bk, 'kernel'); continue; end
                    if strcmpi(bk.name, 'zeros') || ...
                       (isprop(bk,'iszero') && bk.iszero); continue; end
                    if ~isfield(bk.params, 'zk') || isempty(bk.params.zk); continue; end
                    zk = bk.params.zk;
                    [types, coefs] = decompose_kernel_obj(bk);
                    if isempty(types); continue; end
                end

                delta_total = [];
                for kk = 1:numel(types)
                    typ_k = types{kk};
                    coef_k = coefs(kk);
                    if coef_k == 0; continue; end
                    delta = chnk.kernsplit.helm2d_near_correction( ...
                        chnkrlocal(is), chnkrlocal(it), typ_k, zk);
                    if nnz(delta) == 0; continue; end
                    if isempty(delta_total)
                        delta_total = coef_k * full(delta);
                    else
                        delta_total = delta_total + coef_k * full(delta);
                    end
                end
                if isempty(delta_total); continue; end

                rows = offset(it) + r + (0:npts(it)-1)*ndim;
                cols = offset(is) + c + (0:npts(is)-1)*ndim;
                MAT(rows, cols) = MAT(rows, cols) + delta_total;
            end
        end
    end
end
end


function [types, coefs] = decompose_kernel_spec(typ, coef)
% Convert a layout spec (typ, coef) into a list of scalar (type, coef)
% pairs supported by helm2d_near_correction.
typ = lower(typ);
if any(strcmp(typ, {'c','combined'}))
    if numel(coef) ~= 2
        error('add_level1_cross_wlchs: c kernel needs coefs [c1, c2]');
    end
    types = {'d','s'};
    coefs = [coef(1), coef(2)];
elseif any(strcmp(typ, {'sc','spcombined'}))
    if numel(coef) ~= 2
        error('add_level1_cross_wlchs: sc kernel needs coefs [c_sp, c_s]');
    end
    types = {'sp','s'};
    coefs = [coef(1), coef(2)];
else
    types = {typ};
    coefs = coef;
end
end


function [types, coefs] = decompose_kernel_obj(kern)
% Inspect a kernel object and return list of (type, coef) pairs equivalent
% to its bare-kernel evaluation (relevant for cross-edge near correction).
% Returns empty if kernel type is unsupported.
types = {}; coefs = [];
t = lower(kern.type);
switch t
    case {'s','single'}
        c = wlchs_extract_scalar_local(kern, 's');
        types = {'s'}; coefs = c;
    case {'d','double'}
        c = wlchs_extract_scalar_local(kern, 'd');
        types = {'d'}; coefs = c;
    case {'sp','sprime'}
        c = wlchs_extract_scalar_local(kern, 'sp');
        types = {'sp'}; coefs = c;
    case {'dp','dprime','t'}
        c = wlchs_extract_scalar_local(kern, 'dp');
        types = {'dp'}; coefs = c;
    case {'c','combined'}
        if isfield(kern.params,'coefs') && numel(kern.params.coefs) == 2
            types = {'d','s'};
            coefs = [kern.params.coefs(1), kern.params.coefs(2)];
        end
    case {'sc','spcombined'}
        if isfield(kern.params,'coefs') && numel(kern.params.coefs) == 2
            types = {'sp','s'};
            coefs = [kern.params.coefs(1), kern.params.coefs(2)];
        end
    otherwise
        % unsupported (e.g., custom sum-kernel without recognized type)
end
end


function coef = wlchs_extract_scalar_local(kern, t)
% Probe the kernel eval at a fixed src/tgt pair, divide by the analytic
% bare-kernel value to recover the scalar multiplier (handles 2*dkern etc).
zk = kern.params.zk;
src_probe = struct('r',[0;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
tgt_probe = struct('r',[1;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
val = kern.eval(src_probe, tgt_probe);
switch lower(t)
    case 's',  base = 0.25i * besselh(0, zk);
    case 'd',  base = 0.25i * zk * besselh(1, zk);
    case 'sp', base = kernel('helm','sp',zk).eval(src_probe, tgt_probe);
    case 'dp', base = kernel('helm','dp',zk).eval(src_probe, tgt_probe);
    otherwise, base = 1;
end
coef = val / base;
end
