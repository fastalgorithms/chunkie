function val = lap2d_panel_eval(chnkr, type, dens, targets, dlim, eval_opts)
%CHNK.HELSINGO.LAP2D_PANEL_EVAL  off-curve evaluation of Laplace layer
% potentials with kernel-split close correction (Helsing-Ojala / wLCHS).
%
% Conventions: chunkie's S = -log|r|/(2*pi); n_y outward complex normal.
%
% Inputs:
%   chnkr   - chunker (single connected sequence of panels).
%   type    - 's'/'single', 'd'/'double', 'sp'/'sprime' (target-normal
%             derivative of S; requires target normals via struct input).
%   dens    - chnkr.npt x 1 density samples at chunker nodes.
%   targets - 2 x nt real target locations, OR struct with fields .r
%             (2 x nt) and .n (2 x nt) target normals (required for 'sp').
%   dlim    - (optional, default 1.1) close-target threshold.
%   eval_opts - struct with optional opts.forcefmm passed through to
%               chunkerkerneval for the smooth far-field evaluation.

if nargin < 5 || isempty(dlim); dlim = 1.1; end
if nargin < 6; eval_opts = []; end

ngl = chnkr.k;
nch = chnkr.nch;

% Normalize targets spec
if isstruct(targets)
    targ_r = targets.r;
    if isfield(targets,'n'), targ_n = targets.n; else, targ_n = []; end
else
    targ_r = targets;
    targ_n = [];
end

zsrc_all = (chnkr.r(1,:,:) + 1i*chnkr.r(2,:,:));
zsrc_all = zsrc_all(:).';
nz_all   = (chnkr.n(1,:,:) + 1i*chnkr.n(2,:,:));
nz_all   = nz_all(:).';
[~, wgl] = lege.exps(ngl);
awzp_all = chnkr.wts(:).';
d_all    = (chnkr.d(1,:,:) + 1i*chnkr.d(2,:,:));
d_all    = d_all(:).';
wgl_rep  = repmat(wgl(:).', 1, nch);
wzp_all  = d_all .* wgl_rep;
ztgt = targ_r(1,:).' + 1i*targ_r(2,:).';
dens = dens(:);
if ~isempty(targ_n)
    nzt_all = (targ_n(1,:) + 1i*targ_n(2,:)).';
else
    nzt_all = [];
end

% Smooth far-field via chunkerkerneval (FMM-able; avoids dense alloc).
tlow = lower(type);
needs_tgt_normals = any(strcmp(tlow, {'sp','sprime'}));
if needs_tgt_normals && isempty(targ_n)
    error(['CHNK.HELSINGO.LAP2D_PANEL_EVAL: type ''sp'' requires target ' ...
           'normals; pass targets as struct with .n field']);
end
switch tlow
    case {'s','single'},   ker = kernel('lap','s');
    case {'d','double'},   ker = kernel('lap','d');
    case {'sp','sprime'},  ker = kernel('lap','sp');
    otherwise
        error('CHNK.HELSINGO.LAP2D_PANEL_EVAL: unknown type %s', type);
end
ck_opts = struct('forcesmooth', true);
if isstruct(eval_opts) && isfield(eval_opts,'forcefmm') && eval_opts.forcefmm
    ck_opts.forcefmm = true;
end
if needs_tgt_normals
    val = chunkerkerneval(chnkr, ker, dens, struct('r',targ_r,'n',targ_n), ck_opts);
else
    val = chunkerkerneval(chnkr, ker, dens, targ_r, ck_opts);
end

% Per-panel close correction
[rends,~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));
panlen = zeros(nch,1);
for ip = 1:nch
    panlen(ip) = sum(awzp_all((ip-1)*ngl+(1:ngl)));
end
dlim2 = dlim^2;

for ip = 1:nch
    pidx = (ip-1)*ngl + (1:ngl);
    zsrc_p = zsrc_all(pidx);
    nz_p   = nz_all(pidx);
    wzp_p  = wzp_all(pidx);
    awzp_p = awzp_all(pidx);
    a = endsa(ip); b = endsb(ip);

    diff_p = ztgt - zsrc_p;
    d2_p   = real(diff_p).*real(diff_p) + imag(diff_p).*imag(diff_p);
    near   = min(d2_p, [], 2) < dlim2*panlen(ip)^2;
    if ~any(near); continue; end

    ztgt_near = ztgt(near);
    U = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc_p);

    switch lower(type)
        case {'s','single'}
            LogC = chnk.kernsplit.wlchs_target(a, b, ztgt_near, zsrc_p, ...
                nz_p, wzp_p, awzp_p, U, 0);
            Mclose = chnk.kernsplit.lap2d_close('s', ztgt_near, zsrc_p, ...
                nz_p, wzp_p, LogC, []);
            val(near) = val(near) + Mclose * (awzp_p(:).*dens(pidx));

        case {'d','double'}
            [LogC, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt_near, ...
                zsrc_p, nz_p, wzp_p, awzp_p, U, 0);
            Mclose = chnk.kernsplit.lap2d_close('d', ztgt_near, zsrc_p, ...
                nz_p, wzp_p, LogC, CauC);
            val(near) = val(near) + Mclose * dens(pidx);

        case {'sp','sprime'}
            [LogC, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt_near, ...
                zsrc_p, nz_p, wzp_p, awzp_p, U, 0);
            nzt_near = nzt_all(near);
            Mclose = chnk.kernsplit.lap2d_close('sp', ztgt_near, zsrc_p, ...
                nz_p, wzp_p, LogC, CauC, [], nzt_near);
            val(near) = val(near) + Mclose * dens(pidx);
    end
end
end
