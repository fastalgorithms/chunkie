function dmat = diffmat(chnkr, rect_flag)
%DIFFMAT returns the matrix of arc length differentiation along the 
% chunker object 

	if (nargin < 2)
		rect_flag = false;
	end
	tdiff = reshape(lege.dermat(chnkr.k), [chnkr.k, chnkr.k, 1]);
	ds = reshape(vecnorm(chnkr.d), [chnkr.k, 1, chnkr.nch]);
	dmat = (1./ds) .* tdiff;
	if rect_flag
		[~, ~, u, ~] = lege.exps(chnkr.k);
		[~, ~, ~, v] = lege.exps(chnkr.k - 1);
		tmp = v * u(1:end-1, :);
		dmat = pagemtimes(tmp, dmat);
	end
	dmat = num2cell(dmat, [1 2]);	
	dmat = blkdiag(dmat{:});
end
