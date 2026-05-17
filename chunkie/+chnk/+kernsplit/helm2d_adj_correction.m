function [delta, U_src_cell] = helm2d_adj_correction(chnkr, type, zk, U_src_cell)
%CHNK.HELSINGO.HELM2D_ADJ_CORRECTION  adjacent-panel close-correction
% delta for a 2D Helmholtz layer-potential system matrix.
%
% Returns a sparse matrix `delta` such that
%   sysmat_close = sysmat_smooth + delta
% where sysmat_smooth is a system matrix built by simple Gauss-Legendre
% panel quadrature (no kernel-split), and sysmat_close is a properly
% close-corrected version on adjacent panel pairs.  Self and far entries
% of `delta` are zero; only adjacent pairs are populated.
%
% Method: testhelmos's kernel-split.  Log moments computed in the source
% panel's canonical [-1,1] frame via wfrak_log.  Hypersingular moments
% from wlchs_target with ifself=0 (testhelmos uses the same routine for
% on-curve adjacent HypC, since branch-cut issues only affect the log
% moment).  Combined via helm2d_close.
%
% **Supported kernel types** (kernel-specific coefficients in helm2d_close):
%   - 's'  (single layer)
%   - 'd'  (double layer)
%   - 'sp' (S' = normal derivative of single layer, target side)
%   - 'dp' (hypersingular T = D')
% Other types must add a new case in helm2d_close before this routine
% can handle them.  For unsupported types, errors out.
%
% **Limitations:**
%   - Single chunker only (one connected sequence of panels).  Block
%     kernels and multi-edge chunkgraphs handled by the caller.
%   - The target row nearest the shared joint endpoint of two panels
%     loses 3-4 digits relative to other rows (a known wLCHS limitation
%     for targets canonically very close to the panel boundary).  The
%     average operator action over a smooth density typically still
%     achieves 12+ digits.
%
% Inputs:
%   chnkr - chunker object.
%   type  - 's', 'd', or 'dp'.
%   zk    - wavenumber.
%
% Output:
%   delta - sparse npt x npt matrix.  Add to a smooth-GL system matrix
%           to overwrite (sysmat_smooth_adj + delta = wLCHS_close_eval_adj).

ngl = chnkr.k;
nch = chnkr.nch;
npt = chnkr.npt;

[T_GL, wgl, U_GL_panel] = lege.exps(ngl);
% closed form: P_n(-1) = (-1)^n, P_n(1) = 1
P_left  = (-1).^(0:ngl-1);
P_right = ones(1, ngl);

if nargin < 4 || isempty(U_src_cell)
    U_src_cell = cell(1, nch);
end

[rends, ~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));     % nch x 1 (left ends)
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));     % nch x 1 (right ends)
speed_left  = zeros(nch, 1);
speed_right = zeros(nch, 1);
for ip = 1:nch
    di = chnkr.d(:, :, ip);            % 2 x ngl
    cdi = U_GL_panel * (di.');         % ngl x 2 Legendre coefs
    dleft  = P_left  * cdi;            % 1 x 2
    dright = P_right * cdi;
    speed_left(ip)  = norm(dleft);
    speed_right(ip) = norm(dright);
end

ii = []; jj = []; vv = [];

for ip_src = 1:nch
    src_idx = (ip_src-1)*ngl + (1:ngl);
    zsrc = (chnkr.r(1,:,ip_src) + 1i*chnkr.r(2,:,ip_src));
    nz_s = (chnkr.n(1,:,ip_src) + 1i*chnkr.n(2,:,ip_src));
    d_s  = (chnkr.d(1,:,ip_src) + 1i*chnkr.d(2,:,ip_src));
    wzp  = d_s .* wgl(:).';
    awzp = chnkr.wts(:, ip_src).';
    a = endsa(ip_src);
    b = endsb(ip_src);

    % U_src for HypC (curve-specific Legendre Vandermonde inverse)
    if isempty(U_src_cell{ip_src})
        U_src_cell{ip_src} = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc);
    end
    U_src = U_src_cell{ip_src};

    for side = [1, 2]   % side==1: target panel BEFORE source (superdiag)
                        % side==2: target panel AFTER source (subdiag)
        ip_tgt = chnkr.adj(side, ip_src);
        if ip_tgt <= 0; continue; end

        tgt_idx = (ip_tgt-1)*ngl + (1:ngl);
        ztgt = (chnkr.r(1,:,ip_tgt) + 1i*chnkr.r(2,:,ip_tgt)).';
        nzt  = (chnkr.n(1,:,ip_tgt) + 1i*chnkr.n(2,:,ip_tgt)).';

        % alpha = canonical-parameter length ratio of target/source.
        % Computed from speed-ratio at the joint endpoint where the panels
        % meet (their physical d agrees there, so the ratio of canonical
        % d's gives the canonical-length ratio).
        if side == 1
            % target panel BEFORE source: joint = source's LEFT endpoint
            % = target's RIGHT endpoint
            alpha = speed_right(ip_tgt) / speed_left(ip_src);
            trans = -1 - alpha;
        else
            % target panel AFTER source: joint = source's RIGHT endpoint
            % = target's LEFT endpoint
            alpha = speed_left(ip_tgt) / speed_right(ip_src);
            trans = 1 + alpha;
        end
        mscale = alpha;

        [WfrakL, accept] = chnk.kernsplit.wfrak_log(ngl, trans, mscale);
        if isempty(accept); continue; end

        % LogC delta (canonical): WfrakL/W minus smooth-GL log value
        T_acc = T_GL(accept);
        if side == 1
            Ta = T_GL + 1;
            Tb = (1 - T_acc) * alpha;
            log_smooth = log(abs(Ta.' + Tb));
        else
            Ta = T_GL - 1;
            Tb = (T_acc + 1) * alpha;
            log_smooth = log(abs(Ta.' - Tb));
        end
        LogC_w = WfrakL ./ wgl(:).' - log_smooth;
        LogC_full = zeros(ngl, ngl);
        LogC_full(accept, :) = LogC_w;

        % HypC for 'dp' from wlchs_target (ifself=0). Returned regardless of
        % type so we can pass to helm2d_close uniformly. For 's' it's unused.
        if strcmpi(type, 'dp') || strcmpi(type, 't')
            [~, ~, HypC_full] = chnk.kernsplit.wlchs_target( ...
                a, b, ztgt, zsrc, nz_s, wzp, awzp, U_src, 0);
            CauC_full = [];
        elseif strcmpi(type, 'd') || strcmpi(type, 'double') ...
                || strcmpi(type, 'sp') || strcmpi(type, 'sprime')
            [~, CauC_full] = chnk.kernsplit.wlchs_target( ...
                a, b, ztgt, zsrc, nz_s, wzp, awzp, U_src, 0);
            HypC_full = [];
        else
            CauC_full = [];
            HypC_full = [];
        end

        % delta block (M_close * awzp_or_density depending on type)
        switch lower(type)
            case {'s', 'single'}
                Mc = chnk.kernsplit.helm2d_close('s', zk, ztgt, zsrc, nz_s, ...
                                                wzp, LogC_full, []);
                Mblk = Mc .* awzp;       % bake awzp -> M*density convention
            case {'d', 'double'}
                Mblk = chnk.kernsplit.helm2d_close('d', zk, ztgt, zsrc, nz_s, ...
                                                  wzp, LogC_full, CauC_full);
            case {'sp', 'sprime'}
                % Helsing-Karlsson arXiv:1711.09796 §4.4 (eq. KAC):
                % LogC + CauC. CauC vanishes on smooth Γ but is required
                % near corners / branch points.
                Mblk = chnk.kernsplit.helm2d_close('sp', zk, ztgt, zsrc, nz_s, ...
                                                   wzp, LogC_full, CauC_full, [], ...
                                                   nzt, awzp);
            case {'dp', 'dprime', 't'}
                Mblk = chnk.kernsplit.helm2d_close('dp', zk, ztgt, zsrc, nz_s, ...
                                                   wzp, LogC_full, [], HypC_full, ...
                                                   nzt, awzp);
            otherwise
                error('helm2d_adj_correction: unsupported type %s', type);
        end

        % accumulate sparse triples (Mblk IS the delta to add to smooth)
        [I, J] = ndgrid(tgt_idx, src_idx);
        ii = [ii; I(:)];
        jj = [jj; J(:)];
        vv = [vv; Mblk(:)];
    end
end

delta = sparse(ii, jj, vv, npt, npt);
end
