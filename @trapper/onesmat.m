function onesmat = onesmat(trap)

wts = weights(trap);
wts = wts(:);
temp = ones(size(wts));
onesmat = bsxfun(@times,temp,wts.');

end
