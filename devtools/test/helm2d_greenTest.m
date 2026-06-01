helm2d_greenTest0();


function helm2d_greenTest0()
%HELM2D_GREENTEST
%
% test that gradients are set up right in chnk.helm2d.green

iseed = 8675309;
rng(iseed);

zk = randn() + 1i*randn();
start_eps = 1.0;
src = randn(2,1); trg = randn(2,1); mu = randn(2,1); 
srcn = randn(2,1); srcn = srcn/(norm(srcn));

fcn = @(x) fcn1(zk,src,x);
fcnx = @(x) fcnx1(zk,src,x);
fcny = @(x) fcny1(zk,src,x);

niter = 9;
errsf = gradient_check(fcn,trg,start_eps,niter);
errsfx = gradient_check(fcnx,trg,start_eps,niter);
errsfy = gradient_check(fcny,trg,start_eps,niter);

assert(min(errsf) < 1e-10);
assert(min(errsfx) < 1e-8);
assert(min(errsfy) < 1e-8);



end


function [f,g] = fcn1(zk,src,trg)

[f,g] = chnk.helm2d.green(zk,src,trg);

end

function [f,g] = fcnx1(zk,src,trg)

[~,g1,h1] = chnk.helm2d.green(zk,src,trg);
f = g1(1);
g = h1(1:2);

end

function [f,g] = fcny1(zk,src,trg)

[~,g1,h1] = chnk.helm2d.green(zk,src,trg);
f = g1(2);
g = h1(2:3);

end
