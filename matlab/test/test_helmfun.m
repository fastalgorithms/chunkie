%TEST_HELMFUN
%
% test that gradients are set up right in helmfun

addpath('../src','../../mwrap');

zk = randn() + 1i*randn();
start_eps = 1.0;
src = randn(2,1); trg = randn(2,1); mu = randn(2,1); 
srcn = randn(2,1); srcn = srcn/(norm(srcn));

fcn = @(x) fcn1(zk,src,x);
fcnx = @(x) fcnx1(zk,src,x);
fcny = @(x) fcny1(zk,src,x);

gradientTest(fcn,trg,start_eps);
gradientTest(fcnx,trg,start_eps);
gradientTest(fcny,trg,start_eps);

function [f,g] = fcn1(zk,src,trg)

[f,g] = helmfun(zk,src,trg);

end

function [f,g] = fcnx1(zk,src,trg)

[~,g1,h1] = helmfun(zk,src,trg);
f = g1(1);
g = h1(1:2);

end

function [f,g] = fcny1(zk,src,trg)

[~,g1,h1] = helmfun(zk,src,trg);
f = g1(2);
g = h1(2:3);

end