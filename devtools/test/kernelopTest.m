%KERNELOPTEST verify kernel class operations are correct

% setup points
src = []; 
src.r = [[0;0],[0;1]];
src.n = randn(size(src.r)); src.n = src.n./vecnorm(src.n);

targ = []; 
targ.r = [[1;1],[-1;0]]; 
targ.n = randn(size(targ.r)); targ.n = targ.n./vecnorm(targ.n);

% set primitive kernels
skern = kernel('lap','s'); 
dkern = kernel('helm','d',1);

% set scalar
a = pi*1i;

% interleave kernels
fkern1 =kernel([skern;dkern]);

% multiply and divide by a
fkern2 = a.*fkern1;
fkern3 = fkern1/a;

assert(norm(fkern2.eval(src,targ) - a*fkern1.eval(src,targ))<1e-10)
assert(norm(fkern3.eval(src,targ) - 1/a*fkern1.eval(src,targ))<1e-10)

% negate a kernel
nkern = -skern;
assert(norm(nkern.eval(src,targ) + skern.eval(src,targ))<1e-10)

% add and subtract kernels
ckern1 = skern + dkern;
ckern2 = skern - dkern;
assert(norm(ckern1.eval(src,targ) - (skern.eval(src,targ)+dkern.eval(src,targ)))<1e-10)
assert(norm(ckern2.eval(src,targ) - (skern.eval(src,targ)-dkern.eval(src,targ)))<1e-10)



