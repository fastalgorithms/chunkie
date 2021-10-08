function [xs1,wts1,xs0,wts0,ainterp1,ainterp1kron,ainterps0,ainterps0kron]=setuplogquad(k,opdims)
% obtain quadrature nodes and weights, and interpolation matrices for
% logarithmic and nearly logarithmic singularities

qavail = chnk.quadggq.logavail();
[~,i] = min(abs(qavail-k));
assert(qavail(i) == k,'order %d not found, consider using order %d chunks', ...
    k,qavail(i));
[xs1,wts1,xs0,wts0] = chnk.quadggq.getlogquad(k);

ainterp1 = lege.matrin(k,xs1);
temp = eye(opdims(2));
ainterp1kron = kron(ainterp1,temp);

nquad0 = size(xs0,1);

ainterps0kron = zeros(opdims(2)*nquad0,opdims(2)*k,k);
ainterps0 = zeros(nquad0,k,k);

for j = 1:k
    xs0j = xs0(:,j);
    ainterp0_sm = lege.matrin(k,xs0j);
    ainterps0(:,:,j) = ainterp0_sm;
    ainterps0kron(:,:,j) = kron(ainterp0_sm,temp);
end