function onesmat = onesmat(chnkr)

whts = chunkwhts(chnkr);
whts = whts(:);
temp = ones(size(whts));
onesmat = bsxfun(@times,temp,whts.');

end
