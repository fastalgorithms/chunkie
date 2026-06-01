% A simple test to show data fields being used with a kernel.
% The hilbert kernel used here requires an arclength parameterization as
% input, which must be input as a data field.
k = 2; 

% Make a non-circular chnkr
chnkr = chunkerfunc(@(t) [5.*(cos(t))'; (sin(t))'],cparams,pref); 
% Make a data field for the Hilbert kernel
chnkr = chnkr.makedatarows(1);
L = sum(sum(chnkr.wts));
data_tmp = arclengthfun(chnkr)/L;
chnkr.data(:) = data_tmp(:)*pi;

H_mat = chunkermat(chnkr, H_kernel())/L; % The L normalization is necessary

f_1 = sin(2*k*chnkr.data(:)); %data goes from 0 to pi, so k must be even
f_2 = H_mat*f_1; % H_mat should change sin to cos for well-resolved frequencies k.

%This should be zero
norm(f_1.^2 + f_2.^2 - 1)

function ck = H_kernel()
	ck = kernel();
	ck.name = "cotan";
	ck.type = "cot";
	ck.eval = @(s, t) cot_func(s, t);
	ck.opdims = [1,1];
	ck.sing = 'pv';

end

function val = cot_func(s, t)
	theta = (t.data(:)) - (s.data(:).');
	val = cot(theta);
end
