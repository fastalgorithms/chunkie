function wts = weights(trap)
%WEIGHTS

  wts = sqrt(sum((trap.d).^2,1));
  wts = wts*trap.h;

end
