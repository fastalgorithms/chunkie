function varargout = green(prof,dim,src,targ)
%CHNK.RADPP.GREEN a radial Green's function and its derivatives from its
% on-axis profiles
%
% Syntax: [val, grad, hess, ...] = chnk.radpp.green(prof,dim,src,targ)
%
% prof(r) gives the on-axis profiles [G | G_x | G_xx G_yy | G_xxx G_xyy |
% G_xxxx G_xxyy G_yyyy], with x along target - source, which are rotated
% to the derivatives in the format of KERNEL.GREEN.
%
% see also KERNEL.RADPP, CHNK.RADPP.FIT, CHNK.RADPP.FITEVAL

nt = size(targ,2);  ns = size(src,2);
nderiv = max(nargout,1) - 1;

d = cell(1,3);
for i = 1:3
    if i <= size(src,1) && i <= size(targ,1)
        d{i} = reshape(targ(i,:).' - src(i,:),[],1);
    else
        d{i} = zeros(nt*ns,1);
    end
end
r = sqrt(d{1}.^2 + d{2}.^2 + d{3}.^2);
p = prof(r);

u = [d{1:dim}]./r;
u(r == 0,:) = repmat(double(1:dim == 1),nnz(r == 0),1);

varargout{1} = reshape(p(:,1),nt,ns);
if nderiv >= 1
    varargout{2} = reshape(p(:,2).*u,nt,ns,dim);
end
for k = 2:nderiv
    tup = nchoosek(1:dim+k-1,k) - (0:k-1);
    ck = zeros(nt*ns,size(tup,1));
    for m = 1:size(tup,1)
        ck(:,m) = tensval(u,p,tup(m,:));
    end
    varargout{k+1} = reshape(ck,nt,ns,size(ck,2));
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
        error('CHNK.RADPP.GREEN: derivatives of order %d not supported.',k);
end

end
