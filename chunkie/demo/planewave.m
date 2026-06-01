function [y,grad,hess,third] = planewave(kvec,r)

y=exp(1i*sum(bsxfun(@times,kvec(:),r(:,:))));

gx = y*1i*kvec(1);
gy = y*1i*kvec(2);
grad = [gx; gy].';

hxx = y*(1i*kvec(1))^2; 
hxy = y*(1i*kvec(1))*(1i*kvec(2)); 
hyy = y*(1i*kvec(2))^2;
hess = [hxx; hxy; hyy].';

txxx = y*(1i*kvec(1))^3;
txxy = y*(1i*kvec(1))^2*(1i*kvec(2)); 
txyy = y*(1i*kvec(1))*(1i*kvec(2))^2; 
tyyy = y*(1i*kvec(2))^3; 
third = [txxx; txxy; txyy; tyyy].';

end