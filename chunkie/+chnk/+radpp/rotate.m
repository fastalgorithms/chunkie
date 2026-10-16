function comps = rotate(dx,dy,dz,prof,nderiv,dim)
%CHNK.RADPP.ROTATE derivative tensors of a radial function from its
% on-axis profiles
%
% Syntax: comps = chnk.radpp.rotate(dx,dy,dz,prof,nderiv,dim)
%
% dx, dy, dz are displacements target - source (any shape, dz may be 0).
% prof is (numel(dx),nprof), the profiles at r = |(dx,dy,dz)|
%
%    [G | G_x | G_xx G_yy | G_xxx G_xyy | G_xxxx G_xxyy G_yyyy]
%
% x is the direction from the source to the target and y is perpendicular
% to it. The columns are the target derivatives of G on the x-axis, by
% order. Since G is radial, the derivatives with an odd number of y's
% (G_xy, G_xxy, G_yyy, ...) vanish there, and these are the rest. The
% full tensors are obtained from them by rotation.
%
% up to order nderiv <= 4, i.e. nprof = 1, 2, 4, 6, 9.
%
% comps is (numel(dx),ncomp) with the derivatives with respect to the
% target, ordered by derivative order and then as in CHNK.RADPP.TUPLES
%
%   dim = 2 (default), ncomp = 1, 3, 6, 10, 15
%      [G | G_x G_y | H_xx H_xy H_yy | T_xxx T_xxy T_xyy T_yyy | ...]
%   dim = 3, ncomp = 1, 4, 10, 20, 35
%      [G | G_x G_y G_z | H_xx H_xy H_xz H_yy H_yz H_zz | ...]
%
% In both cases r includes dz, so dim = 2 gives the in-plane derivatives
% of a function of the distance in three dimensions.
%
% see also CHNK.RADPP.FIT, CHNK.RADPP.FITEVAL, CHNK.RADPP.TUPLES

if nargin < 6 || isempty(dim), dim = 2; end

dx = dx(:);  dy = dy(:);  dz = dz(:) + 0*dx;
r = sqrt(dx.^2 + dy.^2 + dz.^2);

u = [dx dy dz]./r;
u(r == 0,:) = repmat([1 0 0], nnz(r == 0), 1);
u = u(:,1:dim);

comps = prof(:,1);
if nderiv >= 1
    comps = [comps, prof(:,2).*u];
end

for k = 2:nderiv
    tup = chnk.radpp.tuples(dim,k);
    ck = zeros(numel(r),size(tup,1));
    for m = 1:size(tup,1)
        ck(:,m) = tensval(u,prof,tup(m,:));
    end
    comps = [comps, ck];
end

end


function v = tensval(u,prof,idx)
% one component of the derivative tensor of order numel(idx), in terms of
% the unit vector u and the on-axis profiles
k = numel(idx);
d = @(a,b) double(idx(a) == idx(b));
w = @(a) u(:,idx(a));

switch k
    case 2
        % H_ij = (p-q) u_i u_j + q d_ij
        p = prof(:,3);  q = prof(:,4);
        v = (p-q).*w(1).*w(2) + q*d(1,2);
    case 3
        % T_ijk = (a-3b) u_i u_j u_k + b (d_ij u_k + d_ik u_j + d_jk u_i)
        a = prof(:,5);  b = prof(:,6);
        v = (a-3*b).*w(1).*w(2).*w(3) + ...
            b.*(d(1,2)*w(3) + d(1,3)*w(2) + d(2,3)*w(1));
    case 4
        % F_ijkl = A u_i u_j u_k u_l + B sum_6 d_ij u_k u_l
        %          + C (d_ij d_kl + d_ik d_jl + d_il d_jk)
        a = prof(:,7);  b = prof(:,8);  c = prof(:,9);
        C = c/3;  B = b - C;  A = a - 6*B - 3*C;
        v = A.*w(1).*w(2).*w(3).*w(4) + ...
            B.*(d(1,2)*w(3).*w(4) + d(1,3)*w(2).*w(4) + ...
                d(1,4)*w(2).*w(3) + d(2,3)*w(1).*w(4) + ...
                d(2,4)*w(1).*w(3) + d(3,4)*w(1).*w(2)) + ...
            C*(d(1,2)*d(3,4) + d(1,3)*d(2,4) + d(1,4)*d(2,3));
    otherwise
        error('CHNK.RADPP.ROTATE: derivatives of order %d not supported.',k);
end

end
