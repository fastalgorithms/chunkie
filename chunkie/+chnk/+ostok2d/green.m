%%%%%%%  OSCILLATORY STOKES GREEN FUNCTION    %%%%%%%%%
function [val, varargout] = green(k,src,targ)

% Green function corresponding to = −I∆G_BH + ∇ ⊗ ∇G_BH

[~,~, hess, der3] = chnk.bihar2d.green(k,src,targ);
val = [-hess(:,:,3) hess(:,:,2); hess(:,:,2) -hess(:,:,1)];
% grad_x = [-der3(:,:,3) der3(:,:,2); der3(:,:,2) -der3(:,:,1)];
% grad_y = [-der3(:,:,4) der3(:,:,3); der3(:,:,3) -der3(:,:,2)];

