function [u, grad] = planewave(kx,ky,targ)
x = targ(1,:); y = targ(2,:);
u = exp(1i*(kx*x+ky*y);
grad = zeros(size(targ,2),2);
grad(:,1) = 1i*kx*u;
grad(:,2) = 1i*ky*u;