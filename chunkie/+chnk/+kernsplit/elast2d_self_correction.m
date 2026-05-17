function [delta, U_src_cell] = elast2d_self_correction(chnkr, type, lam, mu, U_src_cell)
%CHNK.KERNSPLIT.ELAST2D_SELF_CORRECTION  self-panel close correction for
% 2D linear elasticity layer potentials, decomposed by sub-block.
%
% Chunkie elasticity convention (chnk.elast2d.kern):
%   beta  = (lam+3mu)/(4 pi mu (lam+2mu))
%   gamma = -(lam+mu)/(4 pi mu (lam+2mu))
%   eta   = mu/(2 pi (lam+2mu))
%   zeta  = (lam+mu)/(pi (lam+2mu))
%
% [S]_ij    = (beta log r + gamma/2) delta_ij + gamma r_i r_j / r^2
% [D]_ij    = -eta (r.n_y)/r^2 delta_ij + eta (r x n_y)/r^2 epsilon_ij ...
%             - zeta r_i r_j (r.n_y)/r^4
%             (epsilon_12 = -1, epsilon_21 = +1)
% [Strac]_ij= +eta (r.n_t)/r^2 delta_ij + eta (r x n_t)/r^2 sigma_ij ...
%             + zeta r_i r_j (r.n_t)/r^4
%             (sigma_12 = +1, sigma_21 = -1)
% [Dalt]_ij = -2 eta (r.n_y)/r^2 delta_ij - zeta r_i r_j (r.n_y)/r^4
%
% Singularity structure of each sub-block on a smooth boundary:
%   S xx,yy  :  log + bounded (uses Laplace S wLCHS + analytic diag limit)
%   S xy=yx  :  bounded only (analytic diag limit)
%   D xx,yy  :  bounded only (analytic diag limit, (r.n_y) gives kappa/2)
%   D xy,yx  :  source-normal Cauchy + bounded
%   Strac xx,yy: bounded only
%   Strac xy,yx: target-normal Cauchy + bounded
%   Dalt all  :  bounded only
%
% Chunkie convention:  awzp_j * n_y_complex = -i * wzp_j  (CCW + outward),
% dn/ds = +kappa tau  giving  (r.n_y)/r^2 -> -kappa/2  along the curve,
%                              (r.n_t)/r^2 -> +kappa/2.
%
% Cauchy close corrections (close minus smooth, per kernel * awzp_j).
% Using chunkie tcau = wzp/(y-x)/i and awzp*n_y = -i*wzp (CCW outward):
%   (r.n_y)/r^2 * awzp_smooth   = -Re(tcau)             --> -Re(CauC)
%   (r x n_y)/r^2 * awzp_smooth = -Im(tcau)             --> -Im(CauC)
%   (r.n_t)/r^2 * awzp_smooth   = -Re(nzt.*tcau./nz)    --> -Re(nzt.*CauC./nz)
%   (r x n_t)/r^2 * awzp_smooth = -Im(nzt.*tcau./nz)    --> -Im(nzt.*CauC./nz)
%
% Inputs:
%   chnkr       - chunker
%   type        - sub-block label, one of
%                 's_xx','s_xy','s_yx','s_yy',
%                 'd_xx','d_xy','d_yx','d_yy',
%                 'strac_xx','strac_xy','strac_yx','strac_yy',
%                 'dalt_xx','dalt_xy','dalt_yx','dalt_yy'
%   lam, mu     - Lame parameters
%   U_src_cell  - cached per-panel inverse Vandermonde from prior calls
%
% Output: sparse npt x npt delta matrix following the chunkermat forcewlchs
%   convention --
%     diagonal entries  = FULL close-eval matrix entry (M_smooth diag is
%                         zeroed by the caller),
%     off-diagonal      = (close - smooth) of the matrix entry,
%   so that  M_smooth + delta  = full close-eval block matrix.

if nargin < 5 || isempty(U_src_cell)
    U_src_cell = cell(1, chnkr.nch);
end

ngl = chnkr.k;
nch = chnkr.nch;
N   = chnkr.npt;
d   = chnkr.d(:,:);
d2  = chnkr.d2(:,:);
w   = chnkr.wts(:);
spd = sqrt(d(1,:).^2 + d(2,:).^2).';
tx  = d(1,:).' ./ spd;
ty  = d(2,:).' ./ spd;
kappa = (d(1,:).*d2(2,:) - d(2,:).*d2(1,:)).' ./ spd.^3;

% material coefficients (match chnk.elast2d.kern)
beta  = (lam + 3*mu) / (4*pi*mu*(lam + 2*mu));
gamma = -(lam + mu) / (4*pi*mu*(lam + 2*mu));
eta   = mu / (2*pi*(lam + 2*mu));
zeta  = (lam + mu) / (pi*(lam + 2*mu));

tlow = lower(type);

switch tlow
    % ------------------------------- S -------------------------------
    case 's_xx'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 's', U_src_cell);
        bounded_diag = (gamma/2 + gamma * tx.^2) .* w;
        delta = (-2*pi*beta) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 's_yy'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 's', U_src_cell);
        bounded_diag = (gamma/2 + gamma * ty.^2) .* w;
        delta = (-2*pi*beta) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case {'s_xy', 's_yx'}
        bounded_diag = gamma * tx .* ty .* w;
        delta = sparse(1:N, 1:N, bounded_diag, N, N);

    % ------------------------------- D -------------------------------
    % D(1,1), D(2,2) include -eta*(r.n_y)/r^2 which is -2*pi*eta * Lap_D;
    % Lap_D's wLCHS captures both diagonal limit and off-diagonal close
    % correction.  Bounded zeta r_i^2 (r.n_y)/r^4 piece keeps the analytic
    % diagonal limit only (matching the Stokes dvel pattern).
    case 'd_xx'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 'd', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* tx.^2 .* w;
        delta = (-2*pi*eta) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 'd_yy'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 'd', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* ty.^2 .* w;
        delta = (-2*pi*eta) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 'd_xy'
        % D(1,2) = +eta (r x n_y)/r^2 - zeta r_x r_y (r.n_y)/r^4.
        % Cauchy piece: +eta * Im(CauC) on off-diag + diag.
        % Bounded piece diag: -zeta tx*ty * (-kappa/2) = +(kappa/2) zeta tx*ty.
        [Cau_delta, U_src_cell] = elast2d_cauchy_self_block( ...
            chnkr, 'src_im', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* tx .* ty .* w;
        delta = (+eta) * Cau_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 'd_yx'
        % D(2,1) = -eta (r x n_y)/r^2 - zeta r_x r_y (r.n_y)/r^4.
        [Cau_delta, U_src_cell] = elast2d_cauchy_self_block( ...
            chnkr, 'src_im', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* tx .* ty .* w;
        delta = (-eta) * Cau_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    % ----------------------------- Strac -----------------------------
    % The eta*(r.n_t)/r^2 piece is Laplace S' up to scale; needs the same
    % off-diagonal Cauchy correction (Re(nzt.*CauC./nz)) as Lap S'.  The
    % zeta r_i r_j (r.n_t)/r^4 piece is smooth on a closed curve and only
    % needs the analytic diagonal limit (matching Stokes strac).
    case 'strac_xx'
        % eta*(r.n_t)/r^2 piece via Lap S'-style correction:
        %   eta*(r.n_t)/r^2 = -2*pi*eta * Lap_S' (Lap_S' = -(r.n_t)/(2 pi r^2)).
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 'sp', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* tx.^2 .* w;
        delta = (-2*pi*eta) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 'strac_yy'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 'sp', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* ty.^2 .* w;
        delta = (-2*pi*eta) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 'strac_xy'
        % Strac(1,2) = +eta (r x n_t)/r^2 + zeta r_x r_y (r.n_t)/r^4.
        [Cau_delta, U_src_cell] = elast2d_cauchy_self_block( ...
            chnkr, 'tgt_im', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* tx .* ty .* w;
        delta = (+eta) * Cau_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 'strac_yx'
        % Strac(2,1) = -eta (r x n_t)/r^2 + zeta r_x r_y (r.n_t)/r^4.
        [Cau_delta, U_src_cell] = elast2d_cauchy_self_block( ...
            chnkr, 'tgt_im', U_src_cell);
        bounded_diag = (kappa/2) .* zeta .* tx .* ty .* w;
        delta = (-eta) * Cau_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    % ----------------------------- Dalt ------------------------------
    case 'dalt_xx'
        % Dalt(1,1) = -2 eta (r.n_y)/r^2 - zeta r_x^2 (r.n_y)/r^4; both
        % bounded.  (r.n_y)/r^2 -> -kappa/2 so diag = +(kappa/2)(2 eta + zeta tx^2).
        diag_v = (kappa/2) .* (2*eta + zeta * tx.^2) .* w;
        delta = sparse(1:N, 1:N, diag_v, N, N);

    case 'dalt_yy'
        diag_v = (kappa/2) .* (2*eta + zeta * ty.^2) .* w;
        delta = sparse(1:N, 1:N, diag_v, N, N);

    case {'dalt_xy', 'dalt_yx'}
        diag_v = (kappa/2) .* zeta .* tx .* ty .* w;
        delta = sparse(1:N, 1:N, diag_v, N, N);

    otherwise
        error('CHNK.KERNSPLIT.ELAST2D_SELF_CORRECTION: type %s not supported', type);
end

end


function [Cau, U_src_cell] = elast2d_cauchy_self_block(chnkr, mode, U_src_cell)
% Build sparse npt x npt close-correction matrix for the SELF block of an
% elasticity Cauchy piece.
%   mode = 'src_im' : delta = Im(CauC)   (source-normal cross piece)
%   mode = 'tgt_im' : delta = Im(nzt .* CauC ./ nz)  (target-normal cross)
%
% On-diagonal entries are the full polynomial wLCHS value (CauC was built
% with ifself = 1 so the diagonal of CauC carries the analytic on-curve
% PV).  Off-diagonal entries are the wLCHS delta over smooth GL.

ngl = chnkr.k;
nch = chnkr.nch;
N   = chnkr.npt;

[~, wgl] = lege.exps(ngl);
[rends, ~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));

ii = []; jj = []; vv = [];
for ip = 1:nch
    idx  = (ip-1)*ngl + (1:ngl);
    zsrc = chnkr.r(1,:,ip) + 1i*chnkr.r(2,:,ip);
    nz   = chnkr.n(1,:,ip) + 1i*chnkr.n(2,:,ip);
    dz   = chnkr.d(1,:,ip) + 1i*chnkr.d(2,:,ip);
    awzp = chnkr.wts(:,ip).';
    wzp  = dz .* wgl(:).';
    a = endsa(ip); b = endsb(ip);

    ztgt = zsrc(:); nzt = nz(:);

    if isempty(U_src_cell{ip})
        U_src_cell{ip} = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc);
    end
    U_src = U_src_cell{ip};

    [~, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt, zsrc, nz, ...
        wzp, awzp, U_src, 1);

    switch mode
        case 'src_im'
            Mblk = -imag(CauC);
        case 'tgt_im'
            Mblk = -imag(nzt .* (CauC ./ nz));
        otherwise
            error('elast2d_cauchy_self_block: unknown mode %s', mode);
    end

    [I, J] = ndgrid(idx, idx);
    ii = [ii; I(:)];
    jj = [jj; J(:)];
    vv = [vv; Mblk(:)];
end
Cau = sparse(ii, jj, vv, N, N);
end
