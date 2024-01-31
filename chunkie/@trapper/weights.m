function wts = weights(trap)
%WEIGHTS - note that this routine must only be used in the
% constructor

  wts = sqrt(sum((trap.d).^2,1));
  wts = wts*trap.h;

end
