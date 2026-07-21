function [val, grad, hess, third, fourth] = quasi_dual_sum(rx,ry,zk,kappa,d)
%CHNK.LAP2DQUAS.QUASI_DUAL_SUM evaluate the quasiperiodic dual (Rayleigh-Bloch)
% sum for the Laplace Green's function using the plane-wave expansion.
%
% Computes the quasiperiodic Laplace Green's function using a
% Fourier-series expansion in the periodic direction, which converges
% rapidly at large |y|.
%
% Syntax: [val,grad,hess,third,fourth] = chnk.lap2dquas.quasi_dual_sum(rx,ry,zk,kappa,d)
%
% Input:
%   rx    - x-components of evaluation points (targets relative to source)
%   ry    - y-components of evaluation points; must be nonzero
%   zk    - wavenumber (set to 0 for Laplace)
%   kappa - (nkappa,1) quasiperiodic phase parameters
%   d     - period (scalar)
%
% Output:
%   val   - (nkappa, npt) Green's function values
%   grad  - (nkappa, npt, 2) gradient [d/drx, d/dry]
%   hess  - (nkappa, npt, 3) Hessian [d^2/drx^2, d^2/drxdry, d^2/dry^2]
%   third - (nkappa, npt, 4) third-order derivatives [xxx, xxy, xyy, yyy]
%   fourth- (nkappa, npt, 5) fourth-order derivatives [xxxx, ..., yyyy]
%
% see also CHNK.LAP2DQUAS.LATTICECOEFS

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