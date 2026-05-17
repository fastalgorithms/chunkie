function Mi = myinv(M, sub_ndim)
% Compute the inverse of a 2x2 block matrix
%
%   M = | A B |
%       | C D |
%
% in the (degenerate or near-singular) case where B*C ~ -I and A is
% close to singular.  The Schur trick on the OTHER block gives a
% well-conditioned formula:
%
%   T   = (A - B*D^{-1}*C)^{-1}
%   M11 = T
%   M12 = -T*B*D^{-1}
%   M21 = -D^{-1}*C*T
%   M22 = D^{-1} + D^{-1}*C*T*B*D^{-1}
%
% Per-node interleaving:
%   sub_ndim = 1 (default, helmos style):  components ordered as
%       (a_1, b_1, a_2, b_2, ...) so A = M(1:2:end, 1:2:end), etc.
%   sub_ndim = d (open-arc Stokes mobility): components ordered as
%       (a_1..a_d, b_1..b_d, a_1..a_d, b_1..b_d, ...) so each "block"
%       in the 2x2 Schur picture is itself a d x d sub-matrix per node.

if nargin < 2 || isempty(sub_ndim), sub_ndim = 1; end

np = size(M, 1);
opdim = 2 * sub_ndim;
Npts  = np / opdim;
assert(mod(np, opdim) == 0, 'myinv: matrix size %d not divisible by 2*sub_ndim=%d', np, opdim);

% Build top vs bottom block index sets.
%   sub_ndim=1 -> top = 1:2:end, bot = 2:2:end (matches the original myinv).
block_off = (0:Npts-1) * opdim;
top_idx = reshape(block_off + (1:sub_ndim).',           1, []);
bot_idx = reshape(block_off + (sub_ndim+1:opdim).',     1, []);

A = M(top_idx, top_idx);
B = M(top_idx, bot_idx);
C = M(bot_idx, top_idx);
D = M(bot_idx, bot_idx);

T   = inv(A - B/D*C);
M11 = T;
M12 = -T*B/D;
M21 = -D\C*T;
M22 = inv(D) + D\C*T*B/D;

Mi = zeros(np);
Mi(top_idx, top_idx) = M11;
Mi(top_idx, bot_idx) = M12;
Mi(bot_idx, top_idx) = M21;
Mi(bot_idx, bot_idx) = M22;
end
