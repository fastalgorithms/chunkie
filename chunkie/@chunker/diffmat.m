function dmat = diffmat(chnkr, op_mode)
%DIFFMAT returns the matrix of arc length differentiation along the 
% chunker object 

	if (nargin < 2)
		op_mode = "square";
	end
	tdiff = reshape(lege.dermat(chnkr.k), [chnkr.k, chnkr.k, 1]);
	ds = reshape(vecnorm(chnkr.d), [chnkr.k, 1, chnkr.nch]);
	dmat = (1./ds) .* tdiff;
	if (op_mode ~= "square")
		[~, ~, u, ~] = lege.exps(chnkr.k);
		[~, ~, ~, v] = lege.exps(chnkr.k - 1);
		tmp = v * u(1:end-1, :);
		dmat = pagemtimes(tmp, dmat);
	end
	if (op_mode == "bvp")
		tmp = zeros(chnkr.k, chnkr.k, chnkr.nch);
		tmp(1:end -1, 1:end, 1:end) = dmat;
		dmat = tmp;
	end
	dmat = num2cell(dmat, [1 2]);	
	dmat = blkdiag(dmat{:});
	if (op_mode == "bvp")
		rh_vec = sum(u, 1);
		u(1:2:end, 1:end) = u(1:2:end, 1:end) * (-1);
		lh_vec = sum(u, 1);
		for i = 1:chnkr.nch
			dmat(chnkr.k*i, (chnkr.k*(i - 1) + 1):chnkr.k*i) = rh_vec;
			j = chnkr.adj(2, i);
			dmat(chnkr.k*i, (chnkr.k*(j - 1) + 1):chnkr.k*j) = lh_vec;
		end
	end
end
