function wts = whts(chnkr)


warning('whts is deprecated and will be removed, use weights instead');

  k = chnkr.k;
  nch = chnkr.nch;
  w = chnkr.wstor;
  wts = reshape(sqrt(sum((chnkr.d).^2,1)),k,nch);
  wts = wts.*w(:);

end
