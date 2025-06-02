chunkerclassunitTest0();


function chunkerclassunitTest0()
%chunkerclassunitTest
%
% setting up some tests for methods of the chunker class
%

% adversarial constructor tests
isuccess = 0;
try
    pref = []; pref.k = -1;
    chnkr = chunker(pref);
catch
    isuccess = 1;
end
assert(isuccess);

isuccess = 0;
try
    pref = []; pref.k = 9; [t,w] = lege.exps(8);
    chnkr = chunker(pref,t,w);
catch
    isuccess = 1;
end
assert(isuccess);

isuccess = 0;
try
    cparams = []; 
    pref = []; pref.k=4; pref.nchmax = 100;
    chnkr = chunkerfunc(@(t) starfish(t),cparams,pref);
catch
    isuccess = 1;
end
assert(isuccess);

cparams = []; cparams.chsmall = 1e-14; cparams.ifclosed=0;
cparams.tb = 2*pi-0.01; cparams.nover=1;
pref = []; pref.k=16; pref.nchmax = 10000;
chnkr = chunkerfunc(@(t) starfish(t),cparams,pref);


for j = 1:chnkr.nch
i1 = chnkr.adj(1,j);
i2 = chnkr.adj(2,j);
if (i1 > 0)
assert(chnkr.adj(2,i1) == j)
end
if (i2 > 0)
assert(chnkr.adj(1,i2) == j)
end
end

% test transforms 
v = [1;2];
com1 = chnkr.r(:,:)*chnkr.wts(:)/sum(chnkr.wts(:));
chnkr2 = v + chnkr;
com2 = chnkr2.r(:,:)*chnkr2.wts(:)/sum(chnkr2.wts(:));
assert(norm(v-(com2-com1))/norm(v) < 1e-14)
chnkr2 = chnkr + v; % add on either side
com2 = chnkr2.r(:,:)*chnkr2.wts(:)/sum(chnkr2.wts(:));
assert(norm(v-(com2-com1))/norm(v) < 1e-14)

A = [1 2; 2 3]; 
chnkr2 = A*chnkr;
assert(abs(area(chnkr2) - det(A)*area(chnkr))/abs(area(chnkr)) < 1e-14);
s = 2; 
chnkr2 = s*chnkr;
assert(abs(area(chnkr2) - s^2*area(chnkr))/abs(area(chnkr)) < 1e-14);
chnkr2 = chnkr*s; % scale on either side 
assert(abs(area(chnkr2) - s^2*area(chnkr))/abs(area(chnkr)) < 1e-14);

caught = false;
try
    chnkr2 = chnkr*A; % only transform by matrix on left
catch
    caught = true;
end
assert(caught);

end


