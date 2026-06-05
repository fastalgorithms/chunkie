%%%%%%%  OSCILLATORY STOKES GREEN FUNCTION    %%%%%%%%%
function val = green(k,src,targ)

% Green function corresponding to = −I∆G_BH + ∇ ⊗ ∇G_BH

% [~,~, hess] = chnk.obihar2d.green(k,src,targ);
[~,~, hess] = chnk.flex2d.helmdiffgreen(k,src,targ);
hess = hess/(k*k);
val = [-hess(:,:,3) hess(:,:,2); hess(:,:,2) -hess(:,:,1)];
