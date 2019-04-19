function tau = taus(chnkr)

k = chnkr.k;
nch = chnkr.nch;
d = chnkr.d;

dd = sqrt(sum(abs(d).^2,1));
tau = bsxfun(@rdivide,chnkr.d,dd);

end

