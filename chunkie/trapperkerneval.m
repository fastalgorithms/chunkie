function fints = trapperkerneval(trap,kern,dens,targs,opts)
%TRAPPERINTKERN compute the convolution of the integral kernel with

wts = weights(trap);
srcinfo = []; srcinfo.r = trap.r; srcinfo.d = trap.d; srcinfo.d2 = trap.d2;
targinfo = []; targinfo.r = targs;
mat = kern(srcinfo,targinfo);
fints = mat*diag(wts)*dens;

end
