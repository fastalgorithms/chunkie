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


end


