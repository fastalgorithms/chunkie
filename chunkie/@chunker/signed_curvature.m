function kappa = signed_curvature(chnkr)
%CURVATURE returns the signed curvature at each node 
%
% (x'y'' - y'x'')/(x'^2+y'^2)^(3/2)
%

assert(size(chnkr.r,1) == 2,'signed curvature only defined in 2D');
kappa = (chnkr.d(1,:).*chnkr.d2(2,:)- ...
    chnkr.d(2,:).*chnkr.d2(1,:))./(sqrt(sum(chnkr.d(:,:).^2,1)).^3);
kappa = reshape(kappa,chnkr.k,chnkr.nch);
end