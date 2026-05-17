function r_f_new = Rfstep(MAT, A, r_f_prev, b_circ, Pbc, PWbc, ...
                          starL, circL, starS, circS, use_myinv, myinv_subdim)
%CHNK.RCIP.RFSTEP one level of the singular-RHS RCIP forward recursion.
%
% Implements the vector form of the R_f recursion (Helsing & Karlsson,
% bgkw6 / J. Comput. Phys. 2022, eq. 28) so that the inversion of the
% previous-level R cancels analytically:
%
%   r_{f,i} = P_W^T * M_i^{-1} * y_i
%
%   M_i = [ A^{-1}            K_b(starL,circL) ]
%         [ K_b(circL,starL)  I + K_b(circL,circL) ]
%
%   y_i(starL) = A^{-1} * r_{f,i-1}        (avoided -- see derivation below)
%   y_i(circL) = b_{ib}^circ
%
% Schur block-inverse on the (starL,starL) block plus cancellation gives:
%
%   VA = K_b(circL,starL) * r_{f,i-1}                 % NB: the A^{-1}*A cancels
%   AU = A * K_b(starL,circL)
%   S  = MAT(circL,circL) - K_b(circL,starL)*A*K_b(starL,circL)
%   inner = S \ (b_{ib}^circ - VA)
%   x(starL) = r_{f,i-1} - AU * inner
%   x(circL) = inner
%
%   r_{f,i}(starS) = PWbc' * x(starL)
%   r_{f,i}(circS) = x(circL)
%
% No explicit inversion of A is needed — the user passes r_{f,i-1} directly
% rather than A^{-1}*r_{f,i-1}.
%
% Inputs:
%   MAT      - 3*k*ndim*nedge x ... type-b system matrix at level i (I+K_b)
%   A        - 2*k*ndim*nedge square: previous R, R_{i-1}
%   r_f_prev - 2*k*ndim*nedge vector: previous r_f^star
%   b_circ   - k*ndim*nedge vector: f at the circL nodes of level-i type-b
%   Pbc, PWbc - prolongation matrices from chnk.rcip.Pbcinit (level-c->b),
%               sizes 2*k*ndim*nedge x k*ndim*nedge.  (Pbc unused here but
%               kept in the signature for symmetry with SchurBana.)
%   starL, circL - bad/good index sets into MAT (size 2*k*ndim*nedge and
%                  k*ndim*nedge respectively)
%   starS, circS - bad/good index sets into the output r_f^star (each
%                  k*ndim*nedge)
%   use_myinv - optional flag: use chnk.rcip.myinv for the Schur complement
%               (matches the open-arc case in SchurBana)
%
% Output:
%   r_f_new  - 2*k*ndim*nedge vector: r_{f,i}^star
%
% See also: chnk.rcip.SchurBana, chnk.rcip.Rcompchunk

if nargin < 11 || isempty(use_myinv)
    use_myinv = false;
end
if nargin < 12 || isempty(myinv_subdim)
    myinv_subdim = 1;
end

K_cs_sL = MAT(circL, starL);
K_sL_cs = MAT(starL, circL);

VA = K_cs_sL * r_f_prev;          % cancels A^{-1}: K_b(circL,starL)*r_{f,i-1}
AU = A * K_sL_cs;                  % R_{i-1} * K_b(starL,circL)
Mschur = MAT(circL, circL) - K_cs_sL * A * K_sL_cs;
if use_myinv
    inner = chnk.rcip.myinv(Mschur, myinv_subdim) * (b_circ - VA);
else
    inner = Mschur \ (b_circ - VA);
end

x_starL = r_f_prev - AU * inner;
x_circL = inner;

r_f_new = zeros(numel(starS) + numel(circS), 1);
r_f_new(starS) = PWbc' * x_starL;
r_f_new(circS) = x_circL;

end
