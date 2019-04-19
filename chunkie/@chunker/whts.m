function wchnk = whts(chnkr)

  k = chnkr.k;
  nch = chnkr.nch;
  [~,w] = legeexps(k);
  wchnk = reshape(sqrt(sum((chnkr.d).^2,1)),k,nch);
  wchnk = wchnk.*bsxfun(@times,w(:),(chnkr.h).');

end
