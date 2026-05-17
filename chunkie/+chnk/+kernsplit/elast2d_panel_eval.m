function val = elast2d_panel_eval(chnkr, type, lam, mu, dens, targets, dlim, eval_opts)
%CHNK.KERNSPLIT.ELAST2D_PANEL_EVAL  off-curve evaluation of 2D linear
% elasticity layer potentials with Wu-Barnett-style kernel-split close
% correction.
%
% Decomposes each elasticity sub-block into a sum of Laplace S/D/S' pieces
% (component-wise) and Stokes svel/dvel pieces (vector-valued), each of
% which is evaluated with its own wLCHS panel-eval that handles both smooth
% far-field and polynomial-exact close correction.  Hence accuracy reaches
% near machine precision down to target-to-curve distances much smaller
% than the panel length, matching the Stokes-mobility precedent.
%
% Decompositions (chnk.elast2d.kern conventions, Lap_S = -log|r|/(2pi),
% Lap_D = (r.n_y)/(2pi r^2), Lap_Sp = -(r.n_t)/(2pi r^2),
% K_svel_ij = (r_i r_j/r^2 + log(1/|r|) δ_ij)/(4pi mu),
% K_dvel_ij = r_i r_j (r.n_y)/(pi r^4)):
%
%   S_ij σ_j  = α_S * Lap_S(σ_i) - π ζ * (K_svel σ)_i + (γ/2) * Σ_i
%               where α_S = -1/(λ+2μ),  Σ_i = integral of σ_i over Γ
%
%   D_ij σ_j  = -2π η * Lap_D(σ_i) - π ζ * (K_dvel σ)_i
%               + η * ∫(r×n_y)/r² σ_j ds * ε_ij  (sum over j, ε_12=+1, ε_21=-1)
%
%   Dalt_ij σ_j = -4π η * Lap_D(σ_i) - π ζ * (K_dvel σ)_i
%
%   Strac_ij σ_j = -2π η * Lap_Sp(σ_i) (target-normal)
%                  - η * ∫(r×n_t)/r² σ_j ds * ε_ij  (sum over j, ε_12=+1, ε_21=-1)
%                  + ζ * ∫ r_i r_j (r.n_t)/r^4 σ_j ds  (no Wu-Barnett yet)
%
% Strac is not yet fully Wu-Barnett because chunkie does not provide a
% Stokes-traction wLCHS panel-eval; its ζ piece uses smooth GL and so
% degrades for very-close targets.  S, D, and Dalt are fully Wu-Barnett.
%
% Inputs:
%   chnkr     - chunker (closed boundary or single connected arc)
%   type      - 's' / 'single', 'd' / 'double', 'strac' / 'straction',
%               'dalt'
%   lam, mu   - Lame parameters
%   dens      - 2 x chnkr.npt or length-2*chnkr.npt vector
%   targets   - 2 x nt real point matrix, or struct with fields .r (2 x nt)
%               and .n (2 x nt) target normals (required for 'strac')
%   dlim      - (optional, default 1.1) close-target panel-length threshold
%   eval_opts - (optional) struct, .forcefmm forwards to chunkerkerneval

if nargin < 7 || isempty(dlim); dlim = 1.1; end
if nargin < 8; eval_opts = []; end

ngl = chnkr.k;
nch = chnkr.nch;
npt = chnkr.npt;

% Normalize target spec
if isstruct(targets)
    targ_r = targets.r;
    if isfield(targets,'n'), targ_n = targets.n; else, targ_n = []; end
else
    targ_r = targets;
    targ_n = [];
end
nt = size(targ_r, 2);

% Normalize density to 2 x npt
if isnumeric(dens) && size(dens,1) == 2
    sigma = dens;
elseif numel(dens) == 2*npt
    sigma = reshape(dens(:), 2, []);
else
    error('CHNK.KERNSPLIT.ELAST2D_PANEL_EVAL: dens must be 2-by-npt or length 2*npt');
end

tlow = lower(type);

% Material coefficients
beta  = (lam + 3*mu) / (4*pi*mu*(lam + 2*mu));
gamma = -(lam + mu) / (4*pi*mu*(lam + 2*mu));
eta   = mu / (2*pi*(lam + 2*mu));
zeta  = (lam + mu) / (pi*(lam + 2*mu));
alpha_S = -2*pi*(beta + gamma);            % = -1/(λ+2μ)
alpha_Ksvel = -pi*zeta;                    % = 4πμγ
alpha_D   = -2*pi*eta;                     % coefficient on Lap_D for D xx,yy
alpha_Dalt = -4*pi*eta;                    % coefficient on Lap_D for Dalt
alpha_Sp   = -2*pi*eta;                    % coefficient on Lap_Sp for Strac
alpha_Kdvel = -pi*zeta;                    % coefficient on K_dvel for D, Dalt

val = zeros(2, nt);

switch tlow
    % ============================ S =============================
    case {'s','single'}
        % Lap_S piece (per component)
        for k = 1:2
            uk = chnk.kernsplit.lap2d_panel_eval(chnkr, 's', ...
                sigma(k,:).', targ_r, dlim, eval_opts);
            val(k,:) = val(k,:) + alpha_S * uk(:).';
        end
        % Stokes svel piece (vector-valued)
        usv = chnk.kernsplit.sto2d_panel_eval(chnkr, 's', mu, sigma, ...
            targ_r, dlim, eval_opts);
        val = val + alpha_Ksvel * usv;
        % γ/2 * Σ constant piece
        awzp_all = chnkr.wts(:).';
        Sig = sigma .* awzp_all;
        Sig_total = sum(Sig, 2);
        val(1,:) = val(1,:) + (gamma/2) * Sig_total(1);
        val(2,:) = val(2,:) + (gamma/2) * Sig_total(2);

    % ============================ D =============================
    case {'d','double'}
        for k = 1:2
            uk = chnk.kernsplit.lap2d_panel_eval(chnkr, 'd', ...
                sigma(k,:).', targ_r, dlim, eval_opts);
            val(k,:) = val(k,:) + alpha_D * uk(:).';
        end
        udv = chnk.kernsplit.sto2d_panel_eval(chnkr, 'd', mu, sigma, ...
            targ_r, dlim, eval_opts);
        val = val + alpha_Kdvel * udv;
        % +η ε_ij ∫(r×n_y)/r² σ_j ds
        u_orth = orthogonal_cauchy_panel_eval(chnkr, sigma, targ_r, ...
            dlim, eval_opts, 'src');
        val(1,:) = val(1,:) + ( eta) * u_orth(2,:);
        val(2,:) = val(2,:) + (-eta) * u_orth(1,:);

    % ========================== Dalt ============================
    case 'dalt'
        for k = 1:2
            uk = chnk.kernsplit.lap2d_panel_eval(chnkr, 'd', ...
                sigma(k,:).', targ_r, dlim, eval_opts);
            val(k,:) = val(k,:) + alpha_Dalt * uk(:).';
        end
        udv = chnk.kernsplit.sto2d_panel_eval(chnkr, 'd', mu, sigma, ...
            targ_r, dlim, eval_opts);
        val = val + alpha_Kdvel * udv;

    % ========================== Strac ===========================
    case {'strac','straction'}
        if isempty(targ_n)
            error(['CHNK.KERNSPLIT.ELAST2D_PANEL_EVAL: ''strac'' requires ' ...
                   'target normals (pass targets as struct with .n)']);
        end
        % Lap_Sp piece (per component, target-normal)
        targ_struct = struct('r', targ_r, 'n', targ_n);
        for k = 1:2
            uk = chnk.kernsplit.lap2d_panel_eval(chnkr, 'sp', ...
                sigma(k,:).', targ_struct, dlim, eval_opts);
            val(k,:) = val(k,:) + alpha_Sp * uk(:).';
        end
        % +η ε_ij ∫(r×n_t)/r² σ_j ds, target-normal flavor
        u_orth = orthogonal_cauchy_panel_eval(chnkr, sigma, targ_r, ...
            dlim, eval_opts, 'tgt', targ_n);
        val(1,:) = val(1,:) + ( eta) * u_orth(2,:);
        val(2,:) = val(2,:) + (-eta) * u_orth(1,:);
        % ζ r_i r_j (r.n_t)/r^4 piece: Stokes-strac panel-eval (full Wu-Barnett).
        % ζ r_i r_j (r.n_t)/r^4 = -π ζ * K_stokes_strac (chunkie convention
        % K_stokes_strac = -(1/π) r_i r_j (r.n_t)/r⁴).
        u_strac_stokes = chnk.kernsplit.sto2d_panel_eval(chnkr, 'strac', mu, ...
            sigma, targ_struct, dlim, eval_opts);
        val = val + (-pi*zeta) * u_strac_stokes;

    otherwise
        error('CHNK.KERNSPLIT.ELAST2D_PANEL_EVAL: unknown type %s', type);
end
end


% ------------------------------------------------------------------
function u_orth = orthogonal_cauchy_panel_eval(chnkr, sigma, targ_r, ...
    dlim, eval_opts, mode, targ_n_in)
% Evaluate v_k(x) = ∫ K(x,y) σ_k(y) ds over Γ where
%   K = (r × n_y)/r²    (mode='src')
%   K = (r × n_t)/r²    (mode='tgt')  -- requires target normals
% with wLCHS close correction (Im(CauC) or Im(nzt·CauC./nz)).

if nargin < 7, targ_n_in = []; end

ngl = chnkr.k;
nch = chnkr.nch;
nt  = size(targ_r, 2);

[~, wgl] = lege.exps(ngl);
zsrc_all = (chnkr.r(1,:,:) + 1i*chnkr.r(2,:,:)); zsrc_all = zsrc_all(:).';
nz_all   = (chnkr.n(1,:,:) + 1i*chnkr.n(2,:,:)); nz_all   = nz_all(:).';
awzp_all = chnkr.wts(:).';
d_all    = (chnkr.d(1,:,:) + 1i*chnkr.d(2,:,:)); d_all    = d_all(:).';
wgl_rep  = repmat(wgl(:).', 1, nch);
wzp_all  = d_all .* wgl_rep;

ztgt = targ_r(1,:).' + 1i*targ_r(2,:).';
if strcmp(mode,'tgt')
    if isempty(targ_n_in)
        error('orthogonal_cauchy_panel_eval: tgt mode requires target normals');
    end
    nzt = (targ_n_in(1,:) + 1i*targ_n_in(2,:)).';
end

% Smooth far-field: sum_j K(x,y_j) σ_j awzp_j (vectorized over targets).
% Compute the 1-component contribution K · σ separately for σ_1 and σ_2.
% K = Im(n_y/δ) (src mode) or Im(n_t/δ) (tgt mode), δ = x - y.
% n_y/δ for src mode: complex 1 x npt per target i.
% Vectorized: build nt x npt matrix of K_ij, then matmul with σ_k.
%
% This is the same cost as chunkerkerneval direct smooth, no FMM.  For
% large nt, npt this dominates; optimize later via chunkerkerneval(lap
% hilb kernel) if needed.
diff_full = ztgt - zsrc_all;                % nt x npt
ny_over_d = nz_all ./ diff_full .* (-1);    % n_y / (x-y) = -n_y/(y-x).
% Wait: we want n_y / δ = n_y / (x - y).  Compute directly:
ny_over_d = nz_all ./ diff_full;            % nt x npt complex
% K_src = Im(n_y/δ), K_tgt = Im(n_t/δ).
switch mode
    case 'src'
        K_full = imag(ny_over_d);
    case 'tgt'
        nt_over_d = nzt ./ diff_full;       % nt x npt complex
        K_full = imag(nt_over_d);
    otherwise
        error('orthogonal_cauchy_panel_eval: unknown mode %s', mode);
end
K_full = K_full .* awzp_all;                % bake in arclength weights

u_orth = zeros(2, nt);
for k = 1:2
    u_orth(k,:) = (K_full * sigma(k,:).').';
end

% Per-panel close correction (poly - smooth)
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
    [~, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt_near, ...
        zsrc_p, nz_p, wzp_p, awzp_p, U, 0);

    % close - smooth correction matrix per source for the integral
    % ∫ K(x,y) σ ds where K is (r × n)/r² (src or tgt normal).
    % Recall: (r × n_y)/r² · awzp_smooth = -Im(tcau); close - smooth = -Im(CauC).
    %         (r × n_t)/r² · awzp_smooth = -Im(nzt·tcau/nz); correction = -Im(nzt·CauC/nz).
    switch mode
        case 'src'
            Mc = -imag(CauC);                       % nt_near x ngl
        case 'tgt'
            Mc = -imag(nzt(near) .* (CauC ./ nz_p));% nt_near x ngl
    end

    sigma_p_1 = sigma(1, pidx).';
    sigma_p_2 = sigma(2, pidx).';
    u_orth(1, near) = u_orth(1, near) + (Mc * sigma_p_1).';
    u_orth(2, near) = u_orth(2, near) + (Mc * sigma_p_2).';
end
end


