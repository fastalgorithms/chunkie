function [delta, U_src_cell] = elast2d_adj_correction(chnkr, type, lam, mu, U_src_cell)
%CHNK.KERNSPLIT.ELAST2D_ADJ_CORRECTION  adjacent-panel close correction
% delta for 2D linear elasticity layer potentials, decomposed by
% sub-block.  See ELAST2D_SELF_CORRECTION for the convention.
%
% Only sub-blocks with a true singular piece (log or Cauchy) need adj
% corrections; bounded pieces are resolved exponentially by smooth GL on
% panel-order-16 adjacent panels and contribute no adj delta.

if nargin < 5 || isempty(U_src_cell)
    U_src_cell = cell(1, chnkr.nch);
end

N = chnkr.npt;

% material coefficients
beta = (lam + 3*mu) / (4*pi*mu*(lam + 2*mu));
eta  = mu / (2*pi*(lam + 2*mu));

tlow = lower(type);

switch tlow
    case {'s_xx', 's_yy'}
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_adj_correction( ...
            chnkr, 's', U_src_cell);
        delta = (-2*pi*beta) * lap_delta;

    case {'s_xy', 's_yx'}
        delta = sparse(N, N);

    case 'd_xx'
        % -eta*(r.n_y)/r^2 piece needs Lap D adj correction; bounded zeta
        % piece resolved by smooth GL.
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_adj_correction( ...
            chnkr, 'd', U_src_cell);
        delta = (-2*pi*eta) * lap_delta;

    case 'd_yy'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_adj_correction( ...
            chnkr, 'd', U_src_cell);
        delta = (-2*pi*eta) * lap_delta;

    case 'd_xy'
        [Cau_adj, U_src_cell] = elast2d_cauchy_adj_block( ...
            chnkr, 'src_im', U_src_cell);
        delta = (+eta) * Cau_adj;

    case 'd_yx'
        [Cau_adj, U_src_cell] = elast2d_cauchy_adj_block( ...
            chnkr, 'src_im', U_src_cell);
        delta = (-eta) * Cau_adj;

    case 'strac_xx'
        % eta*(r.n_t)/r^2 piece via Lap S' adj correction.
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_adj_correction( ...
            chnkr, 'sp', U_src_cell);
        delta = (-2*pi*eta) * lap_delta;

    case 'strac_yy'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_adj_correction( ...
            chnkr, 'sp', U_src_cell);
        delta = (-2*pi*eta) * lap_delta;

    case 'strac_xy'
        [Cau_adj, U_src_cell] = elast2d_cauchy_adj_block( ...
            chnkr, 'tgt_im', U_src_cell);
        delta = (+eta) * Cau_adj;

    case 'strac_yx'
        [Cau_adj, U_src_cell] = elast2d_cauchy_adj_block( ...
            chnkr, 'tgt_im', U_src_cell);
        delta = (-eta) * Cau_adj;

    case {'dalt_xx', 'dalt_yy', 'dalt_xy', 'dalt_yx'}
        delta = sparse(N, N);

    otherwise
        error('CHNK.KERNSPLIT.ELAST2D_ADJ_CORRECTION: type %s not supported', type);
end

end


function [Cau, U_src_cell] = elast2d_cauchy_adj_block(chnkr, mode, U_src_cell)
% Adjacent-panel close-correction for elasticity Cauchy pieces.
%   mode = 'src_im' :  delta = Im(CauC)
%   mode = 'tgt_im' :  delta = Im(nzt .* CauC ./ nz)
%
% Mirrors LAP2D_ADJ_CORRECTION's source/target loop and per-panel-pair
% geometry setup, but extracts Im rather than -Re of CauC (and, for
% target-normal, multiplies by nzt/nz before taking Im).

ngl = chnkr.k;
nch = chnkr.nch;
N   = chnkr.npt;

[T_GL, wgl, U_GL_panel] = lege.exps(ngl);
P_left  = (-1).^(0:ngl-1);
P_right = ones(1, ngl);

[rends, ~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));

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

        [~, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt, zsrc, nz_s, ...
            wzp, awzp, U_src, 0);

        switch mode
            case 'src_im'
                Mblk = -imag(CauC);
            case 'tgt_im'
                Mblk = -imag(nzt .* (CauC ./ nz_s));
            otherwise
                error('elast2d_cauchy_adj_block: unknown mode %s', mode);
        end

        [I, J] = ndgrid(tgt_idx, src_idx);
        ii = [ii; I(:)];
        jj = [jj; J(:)];
        vv = [vv; Mblk(:)];
    end
end
Cau = sparse(ii, jj, vv, N, N);
end
