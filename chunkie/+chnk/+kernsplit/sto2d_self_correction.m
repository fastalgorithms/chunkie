function [delta, U_src_cell] = sto2d_self_correction(chnkr, type, mu, U_src_cell)
%CHNK.KERNSPLIT.STO2D_SELF_CORRECTION  self-panel close correction for the
% 2D Stokes velocity layer-potential kernels, decomposed by sub-block.
%
% Velocity SLP (chunkie convention):
%   K_xx_svel = (1/(4*pi*mu)) (rx^2/r^2 + log(1/r))
%   K_yy_svel = (1/(4*pi*mu)) (ry^2/r^2 + log(1/r))
%   K_xy_svel = K_yx_svel = (1/(4*pi*mu)) rx*ry/r^2
% i.e. a bounded smooth piece (r_i r_j / r^2) plus a Laplace-SLP-times-
% delta_ij log piece.  The log piece equals (1/(2*mu)) times the chunkie
% Laplace SLP kernel; we reuse chnk.kernsplit.lap2d_self_correction for
% that.  The bounded piece needs only the analytic diagonal limit
% tau_i tau_j to finite-fix the smooth-GL diagonal entry.
%
% Velocity DLP (chunkie convention):
%   K_ij_dvel = (1/pi) r_i r_j (n_s . r) / r^4
% bounded everywhere on a smooth closed boundary, with diagonal limit
% kappa tau_i tau_j / (2*pi).  No log; smooth GL on off-diagonals is
% already exponentially accurate at panel order 16, so only the analytic
% diagonal limit is needed.
%
% Inputs:
%   chnkr       - chunker
%   type        - SLP: 'svel_xx', 'svel_xy', 'svel_yx', 'svel_yy'
%                 DLP: 'dvel_xx', 'dvel_xy', 'dvel_yx', 'dvel_yy'
%   mu          - viscosity (used by SLP only; DLP is mu-independent)
%   U_src_cell  - cached per-panel inverse Vandermonde from prior calls
%
% Output: sparse npt x npt delta matrix following the chunkermat
%   forcewlchs convention --
%     diagonal entries  = FULL close-eval matrix entry (M_smooth diag is
%                         zeroed by the caller, so this *replaces* it),
%     off-diagonal      = (close - smooth) of the matrix entry,
% so that  M_smooth + delta  = full close-eval block matrix.

if nargin < 4 || isempty(U_src_cell)
    U_src_cell = cell(1, chnkr.nch);
end
N   = chnkr.npt;
d   = chnkr.d(:,:);
d2  = chnkr.d2(:,:);
w   = chnkr.wts(:);
spd = sqrt(d(1,:).^2 + d(2,:).^2).';
tx  = d(1,:).' ./ spd;
ty  = d(2,:).' ./ spd;
kappa = ((d(1,:).*d2(2,:) - d(2,:).*d2(1,:)).' ./ spd.^3);

switch lower(type)
    case 'svel_xx'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 's', U_src_cell);
        bounded_diag = (1/(4*pi*mu)) * tx.^2 .* w;
        delta = (1/(2*mu)) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case 'svel_yy'
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction( ...
            chnkr, 's', U_src_cell);
        bounded_diag = (1/(4*pi*mu)) * ty.^2 .* w;
        delta = (1/(2*mu)) * lap_delta + sparse(1:N, 1:N, bounded_diag, N, N);

    case {'svel_xy', 'svel_yx'}
        bounded_diag = (1/(4*pi*mu)) * tx.*ty .* w;
        delta = sparse(1:N, 1:N, bounded_diag, N, N);

    case {'dvel_xx', 'strac_xx'}
        % Stokes DLP / traction-of-SLP on a smooth boundary are bounded
        % everywhere with the same diagonal limit -kappa tau_i tau_j /(2*pi).
        %   K_ij_dvel  = +(1/pi) r_i r_j (n_s . r) / r^4   (source normal)
        %   K_ij_strac = -(1/pi) r_i r_j (n_t . r) / r^4   (target normal)
        % Using chunkie's curvature convention dn/ds = +kappa tau, both
        % limits land at  -kappa tau_i tau_j / (2*pi) * w.
        diag_v = -(1/(2*pi)) * kappa .* tx.^2 .* w;
        delta = sparse(1:N, 1:N, diag_v, N, N);

    case {'dvel_yy', 'strac_yy'}
        diag_v = -(1/(2*pi)) * kappa .* ty.^2 .* w;
        delta = sparse(1:N, 1:N, diag_v, N, N);

    case {'dvel_xy', 'dvel_yx', 'strac_xy', 'strac_yx'}
        diag_v = -(1/(2*pi)) * kappa .* tx .* ty .* w;
        delta = sparse(1:N, 1:N, diag_v, N, N);

    otherwise
        error('CHNK.KERNSPLIT.STO2D_SELF_CORRECTION: type %s not supported', type);
end
end
