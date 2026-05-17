function val = helm2d_panel_eval(chnkr, type, zk, dens, targets, dlim, eval_opts)
%CHNK.HELSINGO.HELM2D_PANEL_EVAL  Helmholtz S/D/S' field evaluation
% against a chunker via smooth Gauss-Legendre quadrature with per-panel
% kernel-split close correction (wlchs_target).
%
% Designed to deliver 10+ digits of accuracy at any off-curve target,
% including sub-deepest-panel distances near a corner -- the regime
% where pure adaptive quadrature fails.
%
% Inputs:
%   chnkr   - chunker (single connected sequence of panels). Coordinates
%             must be in the SAME frame as targets (i.e., world coords).
%   type    - 's' / 'single' for Helmholtz single-layer S = (i/4)H_0(zk r);
%             'd' / 'double' for the double-layer D;
%             'sp' / 'sprime' for the target-normal derivative S'
%             (= ∂_{n_t} S) -- requires target normals.
%   zk      - wavenumber.
%   dens    - chnkr.npt x 1 (or row) density samples at chunker nodes.
%   targets - either a 2 x nt real array (works for 's','d') or a struct
%             with fields .r (2 x nt) and .n (2 x nt target normals,
%             required for 'sp'/'sprime').
%   dlim    - (optional, default 1.1) close-target threshold: a target
%             is treated as "close" to a panel if its squared distance
%             to the nearest panel node is < dlim^2 * panel_length^2.
%
% Output:
%   val     - nt x 1 layer potential values at each target.

if nargin < 6 || isempty(dlim); dlim = 1.1; end
if nargin < 7; eval_opts = []; end

% Normalize targets to a 2 x nt array of points; pull out target normals
% if provided (needed for 'sp'/'sprime').
if isstruct(targets)
    if ~isfield(targets,'r')
        error('CHNK.HELSINGO.HELM2D_PANEL_EVAL: targets struct lacks .r');
    end
    tgt_n = [];
    if isfield(targets,'n') && ~isempty(targets.n)
        tgt_n = targets.n;
    end
    targets = targets.r;
else
    tgt_n = [];
end
needs_tgt_normals = any(strcmpi(type, {'sp','sprime'}));
if needs_tgt_normals && isempty(tgt_n)
    error(['CHNK.HELSINGO.HELM2D_PANEL_EVAL: type ''%s'' requires ' ...
           'target normals; pass targets as struct with .n field'], type);
end

ngl = chnkr.k;
nch = chnkr.nch;

% complex coordinates
zsrc_all = (chnkr.r(1,:,:) + 1i*chnkr.r(2,:,:));   % 1 x ngl x nch
zsrc_all = zsrc_all(:).';                           % 1 x npt
nz_all   = (chnkr.n(1,:,:) + 1i*chnkr.n(2,:,:));
nz_all   = nz_all(:).';                             % 1 x npt
% chunker stores d in dr/dt_chunker, where t_chunker is in [-1,1]; the
% Gauss-Legendre weights w_GL are stored as wts/|d|. Here we use the
% standard Legendre weights from lege.exps to reconstruct wzp = z'*w_GL.
[~, wgl] = lege.exps(ngl);
% awzp(i) = |d_i| * wgl_i (this is exactly chnkr.wts) and
% wzp(i)  = d_i (complex) * wgl_i.
awzp_all = chnkr.wts(:).';
d_all    = (chnkr.d(1,:,:) + 1i*chnkr.d(2,:,:));   % 1 x ngl x nch
d_all    = d_all(:).';                              % 1 x npt
wgl_rep  = repmat(wgl(:).', 1, nch);                % 1 x npt
wzp_all  = d_all .* wgl_rep;

ztgt = targets(1,:).' + 1i*targets(2,:).';          % nt x 1

dens = dens(:);

% ---- smooth quadrature for all targets via chunkerkerneval
% (FMM-acceleratable; avoids dense nt x npt allocation at large nt).
switch lower(type)
    case {'s','single'};       ker = kernel('helm','s', zk);
    case {'d','double'};       ker = kernel('helm','d', zk);
    case {'sp','sprime'};      ker = kernel('helm','sprime', zk);
    otherwise
        error('CHNK.HELSINGO.HELM2D_PANEL_EVAL: unknown type %s', type);
end
ck_opts = struct('forcesmooth', true);
if isstruct(eval_opts) && isfield(eval_opts,'forcefmm') && eval_opts.forcefmm
    ck_opts.forcefmm = true;
else
    % accel auto-switches to FMM at nt>200; force dense direct here.
    ck_opts.accel = false;
end
if needs_tgt_normals
    targobj_for_smooth = struct('r', targets, 'n', tgt_n);
else
    targobj_for_smooth = targets;
end
val = chunkerkerneval(chnkr, ker, dens, targobj_for_smooth, ck_opts);

% complex target normals (for 'sp')
if needs_tgt_normals
    nzt_all = (tgt_n(1,:).' + 1i*tgt_n(2,:).');   % nt x 1
end

% ---- per-panel close correction
% panel endpoints in complex form
[rends,~] = chunkends(chnkr);          % dim x 2 x nch
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));     % nch x 1 (left ends)
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));     % nch x 1 (right ends)

panlen = zeros(nch,1);
for ip = 1:nch
    panlen(ip) = sum(awzp_all((ip-1)*ngl+(1:ngl)));
end
panlen2 = (panlen).^2;
dlim2   = dlim^2;

for ip = 1:nch
    pidx = (ip-1)*ngl + (1:ngl);
    zsrc_p = zsrc_all(pidx);
    nz_p   = nz_all(pidx);
    wzp_p  = wzp_all(pidx);
    awzp_p = awzp_all(pidx);
    a      = endsa(ip);
    b      = endsb(ip);

    % screen targets close to this panel
    diff_p = ztgt - zsrc_p;
    d2_p   = real(diff_p).*real(diff_p) + imag(diff_p).*imag(diff_p);
    dnear  = min(d2_p, [], 2);                       % nt x 1
    near   = dnear < dlim2*panlen2(ip);
    if ~any(near); continue; end

    ztgt_near = ztgt(near);
    U = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc_p);

    switch lower(type)
        case {'s','single'}
            LogC = chnk.kernsplit.wlchs_target(a,b,ztgt_near,zsrc_p,nz_p,wzp_p,awzp_p,U,0);
            Mclose = chnk.kernsplit.helm2d_close('s',zk,ztgt_near,zsrc_p,nz_p,wzp_p,LogC,[]);
            val(near) = val(near) + Mclose * (awzp_p(:).*dens(pidx));
        case {'d','double'}
            [LogC, CauC] = chnk.kernsplit.wlchs_target(a,b,ztgt_near,zsrc_p,nz_p,wzp_p,awzp_p,U,0);
            Mclose = chnk.kernsplit.helm2d_close('d',zk,ztgt_near,zsrc_p,nz_p,wzp_p,LogC,CauC);
            val(near) = val(near) + Mclose * dens(pidx);
        case {'sp','sprime'}
            [LogC, CauC] = chnk.kernsplit.wlchs_target(a,b,ztgt_near,zsrc_p,nz_p,wzp_p,awzp_p,U,0);
            nzt_near = nzt_all(near);
            Mclose = chnk.kernsplit.helm2d_close('sp',zk,ztgt_near,zsrc_p,nz_p,wzp_p, ...
                LogC, CauC, [], nzt_near, awzp_p);
            val(near) = val(near) + Mclose * dens(pidx);
    end
end
end
