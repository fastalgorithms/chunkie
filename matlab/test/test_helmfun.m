%TEST_HELMFUN
%
% test that gradients are set up right in helmfun



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

% test vorticity kernel

fvecfun = @(x) doub(zk,src,x,srcn,mu);
vortfun = @(x) doubvort(zk,src,x,srcn,mu);

vortcompare(fvecfun,vortfun,trg,start_eps,6);

fvecfun = @(x) sing(zk,src,x,mu);
vortfun = @(x) singvort(zk,src,x,mu);

vortcompare(fvecfun,vortfun,trg,start_eps,6);

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

function fvec = doub(zk,src,trg,srcn,mu)
mat = ostokeskern(zk,src,trg,srcn,zeros(2,1),'double');
fvec = mat*mu;
end

function vort = doubvort(zk,src,trg,srcn,mu)
cs = 0.0; cd = 1.0;
mat = ostokesvortkern(zk,cs,cd,src,trg,srcn);
vort = mat*mu;
end

function fvec = sing(zk,src,trg,mu)
mat = ostokeskern(zk,src,trg,zeros(2,1),zeros(2,1),'single');
fvec = mat*mu;
end

function vort = singvort(zk,src,trg,mu)
cs = 1.0; cd = 0.0;
mat = ostokesvortkern(zk,cs,cd,src,trg,zeros(2,1));
vort = mat*mu;
end

function vortcompare(fvecfun,vortfun,trg0,hstart,nh)

fvec0 = fvecfun(trg0);
vort0 = vortfun(trg0);
for i = 1:nh
    h = hstart/(10^(i-1));
    trgx = trg0 + [h;0];
    trgy = trg0 + [0;h];
    fvecx = fvecfun(trgx);
    fvecy = fvecfun(trgy);
    
    vortapprox = -(fvecy(1)-fvec0(1))/h + (fvecx(2)-fvec0(2))/h;
    err = abs(vortapprox-vort0);
    fprintf('h %5.2e err %5.2e\n',h,err);
end
end