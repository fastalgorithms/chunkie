datafieldTest0();


function datafieldTest0()
% A simple test to show data fields being used with a kernel.

%% 
% The hilbert kernel used here requires an arclength parameterization as
% input, which must be input as a data field.
k = 5; 

% Make a non-circular chnkr
cparams = []; cparams.nover = 2;
chnkr = chunkerfunc(@(t) starfish(t),cparams);
% Make a data field for the Hilbert kernel
chnkr = chnkr.makedatarows(1);
L = sum(sum(chnkr.wts));
data_tmp = arclengthfun(chnkr)/L;
chnkr.data(:) = data_tmp(:)*2*pi;

opts = []; opt.adaptive_correction = true;
H_mat = chunkermat(chnkr, H_kernel(),opts)/L; % The L normalization is necessary
F = chunkerflam(chnkr,H_kernel(),0);

f_1 = sin(2*k*chnkr.data(:)); %data goes from 0 to pi, so k must be even
f_2 = H_mat*f_1; % H_mat should change sin to cos for well-resolved frequencies k.
f_2_flam = rskelf_mv(F,f_1);

%This should be zero
err1 = norm(f_1.^2 + f_2.^2 - 1)/norm(f_1.^2);

assert(err1 < 1e-8);
assert(norm(f_2-f_2_flam/L)/norm(f_2) < 1e-10);

%% 
% evaluate directional derivative of single layer 
% not an important application but tests some stuff

kern = directional_der_S_kern();
srcinfo = []; srcinfo.r = [5;4];

nt = 10;
targinfo = []; targinfo.r = 0.5*starfish(randn(nt,1)); 
v = [2;1];
targinfo.data = repmat(v,1,nt);

chnkr = chunkerfunc(@(t) starfish(t));

spkern = kernel('l','sp');

unbdry = spkern.eval(srcinfo,chnkr);
sysmat = 0.5*eye(chnkr.npt) + onesmat(chnkr) + chunkermat(chnkr,spkern);

mu = sysmat\unbdry;

deru = chunkerkerneval(chnkr,kern,mu,targinfo);

% try with adaptive integration
opts = []; opts.forceadap = true;
deru2 = chunkerkerneval(chnkr,kern,mu,targinfo,opts);
% try with FLAM
opts = []; opts.flam = true;
deru3 = chunkerkerneval(chnkr,kern,mu,targinfo,opts);

sgradkern = kernel('l','sg');
gradutrue = sgradkern.eval(srcinfo,targinfo);
derutrue = sum(reshape(gradutrue,2,numel(gradutrue)/2).*v,1);

err2 = norm(deru(:)-derutrue(:))/norm(derutrue(:));

assert(err2 < 1e-10);
assert(norm(deru2-deru)/norm(deru) < 1e-10);
assert(norm(deru3-deru)/norm(deru) < 1e-10);


%%



end


function kern = directional_der_S_kern()

kern = kernel();
kern.opdims = [1,1];
kern.sing = 'pv';
kern.eval = @(s,t) directional_der_S(s,t);

end

function submat = directional_der_S(s,t)

v = t.data(:,:);
[~,grad] = chnk.lap2d.green(s.r,t.r);
submat = grad(:,:,1).*(v(1,:).') + grad(:,:,2).*(v(2,:).');

end



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
	val = cot(theta/2);
end
