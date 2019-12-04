function onesmat = onesmat(trap)

whts = whts(trap);
whts = whts(:);
temp = ones(size(whts));
onesmat = bsxfun(@times,temp,whts.');

end
