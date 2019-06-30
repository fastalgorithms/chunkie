function s = arclength(chnkr)

s = bsxfun(@times,squeeze(sum((chnkr.d).^2,1)),chnkr.h(:).');

end