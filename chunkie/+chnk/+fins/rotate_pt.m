function src2 = rotate_pt(src, r, rot)
% rotate src about r by rot
src2 = [];
src2.r = r + rot*(src.r(:,:) - r);

if isfield(src,'n') || isprop(src,'n'),   src2.n = rot*src.r(:,:); end
if isfield(src,'d') || isprop(src,'d'),   src2.d = rot*src.d(:,:); end
if isfield(src,'d2') || isprop(src,'d2'), src2.d2 = rot*src.d2(:,:); end

end