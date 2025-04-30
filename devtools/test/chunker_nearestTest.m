chunker_nearestTest0();


function chunker_nearestTest0()
%chunker.nearest test
%

iseed = 1234;
rng(iseed);

% compare closest computed point on circle geometry to exact answer

circ = chunkerfunc(@(t) chnk.curves.bymode(t,1));
nt = 1000;
thetas = -pi+2*pi*rand(1,nt);  

% the 0.1 here keeps targets away from coordinate singularity of circle
scal = 0.1+2*rand(1,nt);
targs = [cos(thetas); sin(thetas)].*scal;

for j = 1:nt
    [rn,dn,dist,tn,ichn] = nearest(circ,targs(:,j));
    th1 = atan2(targs(2,j),targs(1,j));
    th2 = atan2(rn(2),rn(1));
    err = min(abs(th1-th2),abs(th1-th2+2*pi));
    err = min(err,abs(th1-th2-2*pi));
    assert(err < 1e-12);
end

end


