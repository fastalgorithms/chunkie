function [delta, U_src_cell] = sto2d_adj_correction(chnkr, type, mu, U_src_cell)
%CHNK.KERNSPLIT.STO2D_ADJ_CORRECTION  adjacent-panel close correction
% delta for the 2D Stokes velocity layer potentials, decomposed by
% sub-block.
%
% Stokes svel decomposes as a bounded (r_i r_j / r^2) piece plus a
% (1/(2*mu)) * Laplace_SLP * delta_ij log piece.  On adjacent-panel
% entries the bounded smooth piece is resolved exponentially well by
% smooth Gauss-Legendre quadrature at panel order 16, so the adj
% correction reduces to (1/(2*mu)) * Laplace SLP wLCHS adj correction
% on the (xx, yy) sub-blocks; (xy, yx) need none.
%
% Stokes dvel is fully bounded (no log) on smooth closed curves; smooth
% GL on adjacent panels is exponentially accurate, so all four
% sub-blocks return zero (no adj correction).
%
% Inputs and output convention: see sto2d_self_correction.

if nargin < 4 || isempty(U_src_cell)
    U_src_cell = cell(1, chnkr.nch);
end
N = chnkr.npt;

switch lower(type)
    case {'svel_xx', 'svel_yy'}
        [lap_delta, U_src_cell] = chnk.kernsplit.lap2d_adj_correction( ...
            chnkr, 's', U_src_cell);
        delta = (1/(2*mu)) * lap_delta;

    case {'svel_xy', 'svel_yx'}
        delta = sparse(N, N);

    case {'dvel_xx', 'dvel_xy', 'dvel_yx', 'dvel_yy', ...
          'strac_xx', 'strac_xy', 'strac_yx', 'strac_yy'}
        delta = sparse(N, N);

    otherwise
        error('CHNK.KERNSPLIT.STO2D_ADJ_CORRECTION: type %s not supported', type);
end
end
