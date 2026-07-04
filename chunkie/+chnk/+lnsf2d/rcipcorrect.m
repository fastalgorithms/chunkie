function fcor = rcipcorrect(prm, eta, opts)
%CHNK.LNSF2D.RCIPCORRECT correction handle for the level-wise local
% system matrices inside the RCIP recursion (see the rcip_matcorrect hook
% in chnk.rcip.Rcompchunk).
%
% Returns a function handle
%     MAT = fcor(MAT, chnkrlocal, level, nsub)
% which applies the boundary-layer window quadrature of
% CHNK.LNSF2D.NEARCORRECT to the merged local chunker of the given RCIP
% level, with the kernel scaled by 2 (the identity-normalized convention
% MAT = I + 2(D + eta S) required by the rcip recursion).
%
% Because flagging is based on Im(k) * chunklen, the handle becomes a
% no-op automatically once the recursion has refined the local chunks
% below the boundary-layer scale; deeper levels proceed unmodified.
%
% Notes:
%  * the local chunkers produced by Rcompchunk are recentered at the
%    vertex; the lnsf2d kernels are translation invariant so this is
%    immaterial.
%  * chunkermat on an array of chunkers orders unknowns by concatenation,
%    matching merge(chnkrlocal), so blocks map one-to-one.
%  * cross-edge pairs near the vertex are captured by the distance-based
%    reach criterion in nearcorrect; correcting a pair that did not need
%    it is harmless (the replacement blocks are accurate for any pair).
%
% Syntax:
%   opts.rcip_matcorrect = chnk.lnsf2d.rcipcorrect(prm, eta);
%   sysmat = chunkermat(cgrph, fkern, opts);   % fkern = 2*(D + eta*S)
%
% see also CHNK.LNSF2D.NEARCORRECT, CHNK.RCIP.RCOMPCHUNK

if nargin < 3, opts = []; end
opts.scale = 2.0;

fcor = @(MAT, chnkrlocal, level, nsub) docorrect(MAT, chnkrlocal, ...
    prm, eta, opts);

end

function MAT = docorrect(MAT, chnkrlocal, prm, eta, opts)
chnkrm = merge(chnkrlocal);
MAT = chnk.lnsf2d.nearcorrect(MAT, chnkrm, prm, eta, opts);
end
