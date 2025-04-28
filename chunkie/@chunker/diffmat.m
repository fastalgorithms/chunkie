function dmat = diffmat(chnkr, rect_flag)
	if (nargin < 2)
		rect_flag = false;
	end
	tdiff = reshape(lege.dermat(chnkr.k), [chnkr.k, chnkr.k, 1]);
	ds = reshape(vecnorm(chnkr.d), [chnkr.k, 1, chnkr.nch]);
	dmat = (1./ds) .* tdiff;
	if rect_flag
		[~, ~, u, ~] = lege.exps(chnkr.k);
		[~, ~, ~, v] = lege.exps(chnkr.k - 1);
		tmp = v(:, 1:end-1) * u;
		dmat = pagemtimes(tmp, dmat);
	end
	dmat = num2cell(dmat, [1 2]);	
	dmat = blkdiag(dmat{:})
end
