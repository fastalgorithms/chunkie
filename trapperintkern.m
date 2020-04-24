function fints = trapperintkern(trap,kern,opdims,dens,targs,opts)
%TRAPPERINTKERN compute the convolution of the integral kernel with

wts = weights(trap);
tau = taus(trap);
mat = kern(trap.r,targs,tau,[]);
fints = mat*diag(wts)*dens;

end
