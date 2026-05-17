function [val, grad, hess, third, fourth] = quasi_dual_sum(rx,ry,zk,kappa,d)

kappa = kappa(:).';

tol = 1e-11;
Lbd = sqrt((log(tol))^2/real( min(abs(ry(:))))^2 + real(zk)^2);

rx = rx(:).';
ry = ry(:).';

npt = size(rx,2);

M = ceil(Lbd*d/(2*pi));
ms = reshape((-M:M),1,1,[]);
xi_m = kappa(:) + 2*pi/d*ms;

beta = sqrt(1i*(xi_m-zk)).*sqrt(-1i*(xi_m+zk));

fhat = exp(-beta.*sqrt(ry.^2) + 1i*xi_m.*rx)./(2*beta);
val = sum(fhat,3)/(d);

if nargout > 1
    grad = zeros(length(kappa),npt,2);
    grad(:,:,1) = sum(1i*xi_m.*fhat,3)/d;
    grad(:,:,2) = sum(-beta.*(sqrt(ry.^2)./ry).*fhat,3)/d;
end

if nargout >2
    hess = zeros(length(kappa),npt,3);
    hess(:,:,1) = sum(-xi_m.^2.*fhat,3)/d;
    hess(:,:,2) = sum(-1i*xi_m.*beta.*(sqrt(ry.^2)./ry).*fhat,3)/d;
    hess(:,:,3) = sum((beta.*(sqrt(ry.^2)./ry)).^2.*fhat,3)/d;
end

if nargout > 3
    fac1 = (1i*xi_m);
    fac2 = -beta.*(sqrt(ry.^2)./ry);
    third = zeros(length(kappa),npt,4);
    third(:,:,1) = sum(fac1.^3.*fhat,3)/d;
    third(:,:,2) = sum(fac1.^2.*fac2.*fhat,3)/d;
    third(:,:,3) = sum(fac1.*fac2.^2.*fhat,3)/d;
    third(:,:,4) = sum(fac2.^3.*fhat,3)/d;
end

if nargout > 4
    fac1 = (1i*xi_m);
    fac2 = -beta.*(sqrt(ry.^2)./ry);
    fourth = zeros(length(kappa),npt,5);
    fourth(:,:,1) = sum(fac1.^4.*fhat,3)/d;
    fourth(:,:,2) = sum(fac1.^3.*fac2.*fhat,3)/d;
    fourth(:,:,3) = sum(fac1.^2.*fac2.^2.*fhat,3)/d;
    fourth(:,:,4) = sum(fac1.^1.*fac2.^3.*fhat,3)/d;
    fourth(:,:,5) = sum(fac2.^4.*fhat,3)/d;
end

end