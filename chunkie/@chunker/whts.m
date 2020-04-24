function wchnk = whts(chnkr)


warning('whts is deprecated and will be removed, use weights instead');

  k = chnkr.k;
  nch = chnkr.nch;
  [~,w] = lege.exps(k);
  wchnk = reshape(sqrt(sum((chnkr.d).^2,1)),k,nch);
  wchnk = wchnk.*bsxfun(@times,w(:),(chnkr.h(:)).');

end
