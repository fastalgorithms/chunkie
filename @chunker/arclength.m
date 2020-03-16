function s = arclength(chnkr)

%s = bsxfun(@rdivide,squeeze(sum((chnkr.d).^2,1)),chnkr.h(:).');
s = squeeze(sqrt(sum((chnkr.d).^2,1)));

end