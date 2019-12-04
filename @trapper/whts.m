function wts = whts(trap)

  wts = sqrt(sum((trap.d).^2,1));
  wts = wts*trap.h;

end
