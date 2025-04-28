legeexpsunitTest0();


function legeexpsunitTest0()
%


k = 19;

[x,w,u,v] = lege.exps(k);

dmat = lege.dermat(k,u,v);
imat = lege.intmat(k,u,v);

pv = sin(x);
dpv_true = cos(x);
ipv_true = -cos(x)+cos(-1);

dpv = dmat*pv;
ipv = imat*pv;

assert(norm(dpv-dpv_true) < 1e-12)
assert(norm(ipv-ipv_true) < 1e-14)

cfs = randn(k,1);
cfsint1 = lege.intpol(cfs);
cfsintold = lege.intpol(cfs,'original');
fun = @(t) lege.exev(t,cfs);

for j = 1:k
    itrue = integral(fun,-1,x(j));
    icoefs = lege.exev(x(j),cfsint1);
    assert(abs(itrue-icoefs)< 1e-14);
end



end


