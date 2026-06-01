nv = 3;
z = exp(1j*2*pi*(1:nv)/nv);
verts = [real(z); imag(z)];

opts = [];
opts.lam = 10;
[chnkr,err,err_by_pt] = chnk.smoother.smooth(verts,opts);

assert(err<1E-6)
