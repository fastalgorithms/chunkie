function wts = whts(trap)

warning('whts is deprecated, will be removed. use weights instead');

  wts = sqrt(sum((trap.d).^2,1));
  wts = wts*trap.h;

end
