function [delta, U_src_cell] = lap2d_adj_correction(chnkr, type, U_src_cell)
%CHNK.HELSINGO.LAP2D_ADJ_CORRECTION  adjacent-panel close-correction
% delta for a 2D Laplace layer-potential system matrix.
%
% Returns a sparse npt x npt matrix whose adjacent-panel entries are the
% delta from smooth Gauss-Legendre quadrature.  Self and far entries are
% zero; only adjacent pairs are populated.
%
% Mirrors `helm2d_adj_correction` with the zk -> 0 limit of helm2d_close
% (J_0 -> 1, J_1, J_2 contributions vanish).
%
% Supported types: 's', 'd', 'sp', 'dp'.

ngl = chnkr.k;
nch = chnkr.nch;
npt = chnkr.npt;
tlow = lower(type);
if ~ismember(tlow, {'s','single','d','double','sp','sprime','dp','dprime'})
    delta = sparse(npt, npt);
    if nargout > 1, U_src_cell = {}; end
    return
end

if nargin < 3 || isempty(U_src_cell), U_src_cell = cell(1, nch); end

[T_GL, wgl, U_GL_panel] = lege.exps(ngl);
P_left  = (-1).^(0:ngl-1);
P_right = ones(1, ngl);

[rends, ~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));

% per-panel endpoint speed (for adjacent-panel alpha)
speed_left  = zeros(nch, 1);
speed_right = zeros(nch, 1);
for ip = 1:nch
    di = chnkr.d(:, :, ip);
    cdi = U_GL_panel * (di.');
    dleft  = P_left  * cdi;
    dright = P_right * cdi;
    speed_left(ip)  = norm(dleft);
    speed_right(ip) = norm(dright);
end

ii = []; jj = []; vv = [];

for ip_src = 1:nch
    src_idx = (ip_src-1)*ngl + (1:ngl);
    zsrc = chnkr.r(1,:,ip_src) + 1i*chnkr.r(2,:,ip_src);
    nz_s = chnkr.n(1,:,ip_src) + 1i*chnkr.n(2,:,ip_src);
    d_s  = chnkr.d(1,:,ip_src) + 1i*chnkr.d(2,:,ip_src);
    wzp  = d_s .* wgl(:).';
    awzp = chnkr.wts(:, ip_src).';
    a = endsa(ip_src);
    b = endsb(ip_src);

    if isempty(U_src_cell{ip_src})
        U_src_cell{ip_src} = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc);
    end
    U_src = U_src_cell{ip_src};

    for side = [1, 2]
        ip_tgt = chnkr.adj(side, ip_src);
        if ip_tgt <= 0; continue; end

        tgt_idx = (ip_tgt-1)*ngl + (1:ngl);
        ztgt = (chnkr.r(1,:,ip_tgt) + 1i*chnkr.r(2,:,ip_tgt)).';
        nzt  = (chnkr.n(1,:,ip_tgt) + 1i*chnkr.n(2,:,ip_tgt)).';

        if side == 1
            alpha = speed_right(ip_tgt) / speed_left(ip_src);
            trans = -1 - alpha;
        else
            alpha = speed_left(ip_tgt) / speed_right(ip_src);
            trans = 1 + alpha;
        end
        mscale = alpha;

        [WfrakL, accept] = chnk.kernsplit.wfrak_log(ngl, trans, mscale);
        if isempty(accept); continue; end

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

        % Per-type moment computation
        switch tlow
            case {'s','single'}
                Mc = chnk.kernsplit.lap2d_close('s', ztgt, zsrc, nz_s, ...
                                                wzp, LogC_full);
                Mblk = Mc .* awzp;
            case {'d','double'}
                [~, CauC_full] = chnk.kernsplit.wlchs_target( ...
                    a, b, ztgt, zsrc, nz_s, wzp, awzp, U_src, 0);
                Mblk = chnk.kernsplit.lap2d_close('d', ztgt, zsrc, ...
                    nz_s, wzp, LogC_full, CauC_full);
            case {'sp','sprime'}
                [~, CauC_full] = chnk.kernsplit.wlchs_target( ...
                    a, b, ztgt, zsrc, nz_s, wzp, awzp, U_src, 0);
                Mblk = chnk.kernsplit.lap2d_close('sp', ztgt, zsrc, ...
                    nz_s, wzp, LogC_full, CauC_full, [], nzt);
            case {'dp','dprime'}
                [~, ~, HypC_full] = chnk.kernsplit.wlchs_target( ...
                    a, b, ztgt, zsrc, nz_s, wzp, awzp, U_src, 0);
                Mblk = chnk.kernsplit.lap2d_close('dp', ztgt, zsrc, ...
                    nz_s, wzp, LogC_full, [], HypC_full, nzt);
        end

        [I, J] = ndgrid(tgt_idx, src_idx);
        ii = [ii; I(:)];
        jj = [jj; J(:)];
        vv = [vv; Mblk(:)];
    end
end

delta = sparse(ii, jj, vv, npt, npt);
end
