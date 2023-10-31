function normonesmat = normonesmat(trap)

wts = trap.wts;
rnorms = trap.n;
wts = wts(:);
wts2 = repmat(wts.',2,1);
wts2 = wts2(:).*rnorms(:);

normonesmat = bsxfun(@times,rnorms(:),wts2.');

end
