%TEST_TREFINEMENT
%
% get chunker description of a starfish domain. check if that domain
% satisfies a level restriction in the underlying parameter space. 
% call a refinement routine to fix. 

iseed = 8675309;
rng(iseed);

addpaths_loc();

% define curve

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.25;

% form curve without strict enforcement of level restriction
% in h

start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);
fprintf('%5.2e s : time to build geo\n',t1)
fprintf('originally, nch = %d\n',chnkr.nch)

% check for violations of level restriction

happy = true;
for i = 1:chnkr.nch
    hself = chnkr.h(i);
    i1 = chnkr.adj(1,i);
    i2 = chnkr.adj(2,i);
    h1 = hself; h2 = hself;
    if (i1 > 0) 
        h1 = chnkr.h(i1);
    end
    if (i2 > 0) 
        h2 = chnkr.h(i2);
    end
    if (hself > 2*h1)
        fprintf('oh no! 1 %d\n',i)
        happy = false;
    end
    if (hself > 2*h2)
        fprintf('oh no! 2 %d\n',i)
        happy = false;
    end
end

fprintf('happy? %s\n',mat2str(happy))

% enforce level restriction using refine code with 
% specific options

opts.lvlr = 't'; % refine in t rather than arclength
% standard is 2.05 which is too loose for our purposes
opts.lvlrfac = 1.99; 

start = tic; chnkr = refine(chnkr,opts); t1 = toc(start);
fprintf('%5.2e s : time to refine geo\n',t1)
fprintf('now, nch = %d\n',chnkr.nch)

happy = true;
for i = 1:chnkr.nch
    hself = chnkr.h(i);
    i1 = chnkr.adj(1,i);
    i2 = chnkr.adj(2,i);
    h1 = hself; h2 = hself;
    if (i1 > 0) 
        h1 = chnkr.h(i1);
    end
    if (i2 > 0) 
        h2 = chnkr.h(i2);
    end
    if (hself > 2*h1)
        fprintf('oh no! 1 %d\n',i)
        happy = false;
    end
    if (hself > 2*h2)
        fprintf('oh no! 2 %d\n',i)
        happy = false;
    end
end

fprintf('happy? %s\n',mat2str(happy))
