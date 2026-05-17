function val = sto2d_panel_eval(chnkr, type, mu, dens, targets, dlim, eval_opts)
%CHNK.KERNSPLIT.STO2D_PANEL_EVAL  off-curve evaluation of Stokes velocity
% layer potentials with kernel-split close correction.
%
% Following Wu, Barnett et al. (arXiv 1909.00049) eq (32)-(34):
%   S sigma = (1/(2*mu))[(S sigma_1, S sigma_2)
%                        + grad S(y . sigma)
%                        - x_1 grad S sigma_1
%                        - x_2 grad S sigma_2]
%   D sigma = ((1/(2*pi)) Re(I_C tau_1, I_C tau_2)
%              + grad D(y . sigma)
%              - x_1 grad D sigma_1
%              - x_2 grad D sigma_2)
% where tau_1 = (sigma_1+i*sigma_2)*Re(n_y)/n_y, tau_2 likewise with Im.
% The Laplace S/D potentials use chunkie's G = -log r/(2*pi) convention,
% and grad denotes target gradient.
%
% Inputs:
%   chnkr   - chunker (single connected sequence of panels).
%   type    - 's' / 'single' (Stokes velocity SLP) or 'd' / 'double'.
%   mu      - Stokes viscosity.
%   dens    - 2 * chnkr.npt density samples (interleaved [s1 s2 s1 s2 ...]
%             OR 2-by-npt grid).
%   targets - 2 x nt real target locations.
%   dlim    - (optional, default 1.1) close-target threshold.
%   eval_opts - (optional) struct with .forcefmm passed to chunkerkerneval.
%
% Output:
%   val     - 2 x nt array of Stokes velocity at each target.

if nargin < 6 || isempty(dlim); dlim = 1.1; end
if nargin < 7; eval_opts = []; end

ngl = chnkr.k;
nch = chnkr.nch;
npt = chnkr.npt;

% normalize density to 2 x npt
if size(dens, 1) == 2
    sigma = dens;
elseif numel(dens) == 2*npt
    sigma = reshape(dens(:), 2, []);
else
    error('CHNK.KERNSPLIT.STO2D_PANEL_EVAL: dens must be 2-by-npt or length 2*npt');
end

% complex source / target data
zsrc_all = chnkr.r(1,:,:) + 1i*chnkr.r(2,:,:); zsrc_all = zsrc_all(:).';  % 1 x npt
nz_all   = chnkr.n(1,:,:) + 1i*chnkr.n(2,:,:); nz_all   = nz_all(:).';
[~, wgl] = lege.exps(ngl);
awzp_all = chnkr.wts(:).';
d_all    = chnkr.d(1,:,:) + 1i*chnkr.d(2,:,:); d_all    = d_all(:).';
wgl_rep  = repmat(wgl(:).', 1, nch);
wzp_all  = d_all .* wgl_rep;
y_all    = chnkr.r(:,:);                             % 2 x npt (real)

if isstruct(targets)
    targ_r = targets.r;
    if isfield(targets,'n'), targ_n = targets.n; else, targ_n = []; end
else
    targ_r = targets;
    targ_n = [];
end
ztgt = targ_r(1,:).' + 1i*targ_r(2,:).';              % nt x 1
nt   = size(targ_r, 2);
val  = zeros(2, nt);
if ~isempty(targ_n)
    nzt_all = (targ_n(1,:) + 1i*targ_n(2,:)).';       % nt x 1
else
    nzt_all = [];
end

% Note: we do NOT use chunkie's direct K_stokes kernel for the smooth
% far-field, because the decomposition (Cauchy + Laplace D gradients)
% used for the close correction is only equal to K_stokes after
% INTEGRATION (eq (31) of Wu et al.).  Per-source the two forms differ;
% mixing them produces wrong values at near-singular targets.  Instead,
% the smooth far-field below is built panel-by-panel using the same
% decomposition form, with smooth GL replacing wLCHS for far panels.

if ~strcmpi(type,'s') && ~strcmpi(type,'single') && ...
   ~strcmpi(type,'d') && ~strcmpi(type,'double') && ...
   ~strcmpi(type,'strac') && ~strcmpi(type,'straction')
    error('CHNK.KERNSPLIT.STO2D_PANEL_EVAL: unknown type %s', type);
end

is_strac = strcmpi(type,'strac') || strcmpi(type,'straction');
if is_strac && isempty(nzt_all)
    error(['CHNK.KERNSPLIT.STO2D_PANEL_EVAL: ''strac'' requires target ' ...
           'normals; pass targets as struct with .r and .n']);
end

% --- per-panel evaluation, all panels (smooth + close correction in
% the same decomposition form).
[rends,~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));
panlen = zeros(nch,1);
for ip = 1:nch
    panlen(ip) = sum(awzp_all((ip-1)*ngl+(1:ngl)));
end
dlim2 = dlim^2;

twopi = 2*pi;
xt1 = targ_r(1,:).';
xt2 = targ_r(2,:).';

for ip = 1:nch
    pidx = (ip-1)*ngl + (1:ngl);
    zsrc_p = zsrc_all(pidx);
    nz_p   = nz_all(pidx);
    wzp_p  = wzp_all(pidx);
    awzp_p = awzp_all(pidx);
    y_p    = y_all(:, pidx);
    s1_p   = sigma(1, pidx);
    s2_p   = sigma(2, pidx);
    a = endsa(ip); b = endsb(ip);

    diff_p = ztgt - zsrc_p;
    d2_p   = real(diff_p).*real(diff_p) + imag(diff_p).*imag(diff_p);
    is_near = min(d2_p, [], 2) < dlim2*panlen(ip)^2;

    % I_C and I_H smooth-GL per-source moments (nt x ngl) for THIS panel.
    % Convention: I_C(tau) = integral tau(y)/(y-x) dy.  Smooth GL is
    % sum_j wzp_j/(zsc_j - ztg_i) tau_j = sum_j wzp_j/(-diff_p) tau_j.
    IC_smooth = -wzp_p ./ diff_p;        % nt x ngl  (note minus sign)
    IH_smooth =  wzp_p ./ diff_p.^2;     % I_H = int tau/(y-x)^2 dy, even power

    % For close targets, replace smooth with polynomial-evaluated value
    if any(is_near)
        ztgt_n = ztgt(is_near);
        if strcmpi(type,'s') || strcmpi(type,'single')
            U = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc_p);
            [LogC, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt_n, ...
                zsrc_p, nz_p, wzp_p, awzp_p, U, 0);
        else
            U = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc_p);
            [LogC, CauC, HypC] = chnk.kernsplit.wlchs_target(a, b, ztgt_n, ...
                zsrc_p, nz_p, wzp_p, awzp_p, U, 0);
        end
        % wcmpC = (poly - smooth)/1i for I_C per source.  So:
        %    poly_per_source = i*wcmpC + smooth
        %    delta = poly - smooth = i*wcmpC
        % Add deltas onto the smooth values for near rows:
        IC_smooth(is_near, :) = IC_smooth(is_near, :) + 1i*CauC;
        if strcmpi(type,'d') || strcmpi(type,'double') || is_strac
            % Same for I_H: poly = i*wcmpH + smooth, delta = i*wcmpH
            IH_smooth(is_near, :) = IH_smooth(is_near, :) + 1i*HypC;
        end
        % Also need close-corrected I_L for SLP via Laplace S(tau): the
        % decomposition for SLP needs S(s1), S(s2), S(y.sigma) and grad S
        % of same.  S corrections use LogC; grad S uses I_C/(in_y).
    end

    % --- compute the decomposition contributions per target
    ydots = y_p(1,:).*s1_p + y_p(2,:).*s2_p;     % 1 x ngl, (y.sigma)

    switch lower(type)
        case {'s','single'}
            % S sigma = (1/(2*mu))[(S s1, S s2)
            %                      + grad S(y.sigma) - x_1 grad S(s1) - x_2 grad S(s2)]
            % S(tau) per-source moment (smooth or poly-corrected) is via
            % Laplace SLP = -1/(2pi) I_L.  We use the same wlchs LogC
            % moment as in lap2d (for far targets, the smooth_S_per_src is
            % just -log(r)*awzp/(2pi); for close, add the LogC delta).
            log_r_smooth = log(abs(diff_p));        % nt x ngl
            % Laplace S smooth per source: -log(r)/(2pi)*awzp
            S_per_src = -log_r_smooth .* awzp_p / twopi;
            if any(is_near)
                % Close correction adds -LogC*awzp/(2pi)
                S_per_src(is_near, :) = S_per_src(is_near, :) + ...
                    (-LogC .* awzp_p / twopi);
            end
            % grad S per source: x-component = (1/(2pi)) Re(I_C(tau/(in_y))_per_src)
            % = (1/(2pi)) Re(IC_per_src/(i*nz_j)) when tau is the
            % Lagrange basis at j.  But here tau is generic; per-source j
            % the integrand contribution to I_C(tau/(in_y)) at row i is
            % IC_smooth[i,j] / (i*nz_j) * tau_j.  So grad-S coefficient for
            % source j of value tau_j is:
            %   d/dx1 = (1/(2pi)) Re(IC_smooth[i,j] / (i*nz_j))
            %   d/dx2 = -(1/(2pi)) Im(IC_smooth[i,j] / (i*nz_j))
            G_per_src = IC_smooth ./ (1i * nz_p) / twopi;     % complex, nt x ngl
            GxRe = real(G_per_src);
            GxIm = -imag(G_per_src);

            ux = (S_per_src * s1_p(:)) ...
                + (GxRe * ydots(:)) ...
                - xt1 .* (GxRe * s1_p(:)) ...
                - xt2 .* (GxRe * s2_p(:));
            uy = (S_per_src * s2_p(:)) ...
                + (GxIm * ydots(:)) ...
                - xt1 .* (GxIm * s1_p(:)) ...
                - xt2 .* (GxIm * s2_p(:));
            ux = ux / (2*mu);
            uy = uy / (2*mu);

        case {'d','double'}
            % D sigma decomposition (Wu et al. eq 31, in chunkie's
            % per-source "smooth GL" form).  First term:
            %   (1/(2pi)) integral n_y/rho^2 (r.sigma) ds
            % which is REAL.  In complex form using the I_C moment, this
            % becomes (1/(2pi)) (-Im) of I_C applied to a complex tau.
            % Paper (arXiv:1909.00049) eq 33 defines tau_k = (sigma_1
            % + i sigma_2) * Re(n_y)/n_y for k=1 and Im for k=2; the
            % paper writes "Re" but empirical verification shows the
            % correct extraction is -Im (this is just a sign-convention
            % difference, possibly a typo in the paper).
            sigma_c = s1_p + 1i*s2_p;
            tau1_p  = sigma_c .* (real(nz_p) ./ nz_p);
            tau2_p  = sigma_c .* (imag(nz_p) ./ nz_p);
            ICt1 = IC_smooth * tau1_p(:);
            ICt2 = IC_smooth * tau2_p(:);

            % grad D per source j:
            %   d/dx1 D(tau) = -(1/(2pi)) Im(IH_per_src) tau_j
            %   d/dx2 D(tau) = -(1/(2pi)) Re(IH_per_src) tau_j
            IH_yS = IH_smooth * ydots(:);
            IH_s1 = IH_smooth * s1_p(:);
            IH_s2 = IH_smooth * s2_p(:);

            ux = -imag(ICt1)/twopi ...
                + (-imag(IH_yS)/twopi) ...
                + xt1 .* (imag(IH_s1)/twopi) ...
                + xt2 .* (imag(IH_s2)/twopi);
            uy = -imag(ICt2)/twopi ...
                + (-real(IH_yS)/twopi) ...
                + xt1 .* (real(IH_s1)/twopi) ...
                + xt2 .* (real(IH_s2)/twopi);

        case {'strac','straction'}
            % Stokes traction-of-SLP kernel K_strac = -(1/π) r_i r_j (r·n_t)/r⁴
            % is the same shape as K_dvel = (1/π) r_i r_j (r·n_y)/r⁴ but with
            % target normal n_t (constant per target) replacing source normal
            % n_y, and overall sign flipped.  We obtain the K_strac
            % Wu-Barnett decomposition by substituting in the 'd' code:
            %   tau_k:  use Re(n_t(x))/n_y(y)  and  Im(n_t(x))/n_y(y)
            %   IH[g]:  use g/n_y on source side, multiply by n_t(x) at target
            % then negate the entire result.
            sigma_c = s1_p + 1i*s2_p;
            nzt_near = nzt_all(:);            % nt x 1, all targets
            % sources factored by 1/n_y (chunkie tcau/thyp bake n_y in via
            % awzp*n_y = -i*wzp, so dividing source-side density by n_y
            % cancels that and lets us multiply by n_t at the target).
            sigma_c_over_nz   = sigma_c ./ nz_p;
            ydots_over_nz     = ydots ./ nz_p;
            s1_over_nz        = s1_p ./ nz_p;
            s2_over_nz        = s2_p ./ nz_p;

            ICt0 = IC_smooth * sigma_c_over_nz(:);          % nt x 1 complex
            ICt1 = real(nzt_near) .* ICt0;                  % τ_1_strac · IC
            ICt2 = imag(nzt_near) .* ICt0;

            IH_yS = nzt_near .* (IH_smooth * ydots_over_nz(:));
            IH_s1 = nzt_near .* (IH_smooth * s1_over_nz(:));
            IH_s2 = nzt_near .* (IH_smooth * s2_over_nz(:));

            % K_strac = -K_dvel-with-n_t, so flip every sign of the 'd' result.
            ux = +imag(ICt1)/twopi ...
                + (+imag(IH_yS)/twopi) ...
                - xt1 .* (imag(IH_s1)/twopi) ...
                - xt2 .* (imag(IH_s2)/twopi);
            uy = +imag(ICt2)/twopi ...
                + (+real(IH_yS)/twopi) ...
                - xt1 .* (real(IH_s1)/twopi) ...
                - xt2 .* (real(IH_s2)/twopi);
    end

    val(1, :) = val(1, :) + ux.';
    val(2, :) = val(2, :) + uy.';
end
end
