function kappa = signed_curvature(cgrph)
%CURVATURE returns the signed curvature at each node 
%
% (x'y'' - y'x'')/(x'^2+y'^2)^(3/2)
%

assert(size(cgrph.r,1) == 2,'signed curvature only defined in 2D');
kappa = (cgrph.d(1,:).*cgrph.d2(2,:)- ...
    cgrph.d(2,:).*cgrph.d2(1,:))./(sqrt(sum(cgrph.d(:,:).^2,1)).^3);
nch = sum(horzcat(cgrph.echnks.nch));
kappa = reshape(kappa,cgrph.k,nch);
end