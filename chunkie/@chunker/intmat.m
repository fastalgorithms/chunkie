function imat = intmat(chnkr)

	panel_int = reshape(lege.intmat(chnkr.k), [chnkr.k, chnkr.k, 1]);
	ds = reshape(vecnorm(chnkr.d), [1, chnkr.k, chnkr.nch]);
	panel_int = panel_int .* ds;
	[~, w, ~, ~] = lege.exps(chnkr.k);
	ds = reshape(w, [chnkr.k, 1]) .* reshape(ds, [chnkr.k, chnkr.nch]);
	% Set diagonal blocks
	panel_int = num2cell(panel_int, [1 2]);
	imat = blkdiag(panel_int{:});
	% maybe not the fast way to do this, but find the panel ordering
	order_vec = zeros(chnkr.nch, 1);
	j = 1;
	for i = 1:chnkr.nch
		order_vec(i) = j;
		j = chnkr.adj(2, j);	
	end
	% Set off diagonal matrix components
	for i = 2:chnkr.nch
		for j = 1:(i - 1)
			
			u = 1 + (order_vec(i) - 1)*chnkr.k;
			v = 1 + (order_vec(j) - 1)*chnkr.k;
			tmp = ones(chnkr.k, 1) * ds(1:chnkr.k, order_vec(j))';
			imat(u:u + chnkr.k - 1, v:v + chnkr.k - 1) = tmp;
		end
	end
end
