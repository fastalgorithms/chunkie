
%FLAGNEARTEST
%
% This file tests the flagnear routine against brute force

%

clearvars; close all;
addpaths_loc();
cparams = [];
cparams.eps = 1.0e-6;
pref = []; 
pref.k = 16;
narms = 10;
amp = 0.5;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0,'adjacency issues after chunk build starfish');


fac = 0.7;

nt = 1000;
scal = 2*rand(1,nt);
tr = 2*pi*rand(1,nt);

targs = bsxfun(@times,starfish(tr,narms,amp),scal);

opts = []; opts.fac = fac;
start = tic; flag = flagnear(chnkr,targs,opts); t1 = toc(start);
fprintf('%5.2e s : time for flagnear routine\n',t1);

% brute force (almost intentionally dumb)

lens = chunklen(chnkr);

flag2 = sparse([],[],[],nt,chnkr.nch);
for i = 1:chnkr.nch
    ris = chnkr.r(:,:,i);
    leni = lens(i);
    for j = 1:chnkr.k
        rj = ris(:,j);
        for l = 1:nt
            dist = sqrt(sum((rj-targs(:,l)).^2,1));
            if dist < fac*leni
                flag2(l,i) = true;
            end
        end
    end
end

assert(nnz(flag2 ~= flag) == 0)

    
