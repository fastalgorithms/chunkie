function f = profiles(green,r,nderiv,dim)
%CHNK.RADPP.PROFILES on-axis profiles of a radial Green's function
%
% Syntax: f = chnk.radpp.profiles(green,r,nderiv,dim)
%
% green is called as [val, grad, hess, ...] = green(src,targ) with a
% source at the origin and targets on the x-axis, in dim dimensions. See
% CHNK.RADPP.UNPACK. f is (numel(r),nprof), with columns
%
%    [G | G_x | G_xx G_yy | G_xxx G_xyy | G_xxxx G_xxyy G_yyyy]
%
% x is the direction from the source to the target and y is perpendicular
% to it. The columns are the target derivatives of G on the x-axis, by
% order. Since G is radial, the derivatives with an odd number of y's
% (G_xy, G_xxy, G_yyy, ...) vanish there, and these are the rest. The
% full tensors are obtained from them by rotation.
%
% up to order nderiv <= 4.
%
% see also CHNK.RADPP.ROTATE, CHNK.RADPP.UNPACK

if nargin < 4 || isempty(dim), dim = 2; end

% the components used from the derivatives of each order
sel = {1, [1 1; 2 2], [1 1 1; 1 2 2], [1 1 1 1; 1 1 2 2; 2 2 2 2]};

r = r(:);
n = numel(r);
out = cell(1,nderiv+1);
[out{:}] = green(zeros(dim,1),[r.'; zeros(dim-1,n)]);

f = reshape(out{1},[],1);
for k = 1:nderiv
    c = chnk.radpp.unpack(out{k+1},dim,k);
    tup = chnk.radpp.tuples(dim,k);
    for j = 1:size(sel{k},1)
        m = find(all(tup == sel{k}(j,:),2));
        f = [f, c(:,m)];
    end
end

end
