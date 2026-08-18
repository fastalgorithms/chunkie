function y = mydepth(zsctr,x,Ng)
%CHNK.HELSINGO.MYDEPTH  evaluate a polynomial that interpolates the
% imaginary parts of the recentered source nodes (with imag=0 forced at
% +/-1) at a real query point x.
%
% Used by wlchs_target to detect the convex/concave panel exception
% when deciding whether a target lies on the "left" or "right" of the
% complexified contour.

A  = ones(Ng+2);
xx = [-1; real(zsctr(:)); 1];
for k = 2:Ng+2
    A(:,k) = xx.*A(:,k-1);
end
coeff = A \ [0; imag(zsctr(:)); 0];
y = coeff(Ng+2);
for k = Ng+1:-1:1
    y = y*x + coeff(k);
end
end
