rcipTest0();

function rcipTest0()
%RCIPTEST
%
% This file tests the rcip routines for solving the exterior dirichlet 
% problem on a domain defined by two arcs of circles meeting at two vertices


ncurve = 2;
chnkr(1,ncurve) = chunker();

% set wave number
zk = 1.1;

funs = cell(2,1); 
funs{1} = @(t) circle(t); funs{2} = @(t) circle(t);
cparamcell = cell(2,1);
cparams1 = []; cparams1.ta = 0; cparams1.tb = pi/1.1; cparams1.nchmin = 1;
cparams2 = []; cparams2.ta = 0; cparams2.tb = pi/1.4; cparams2.nchmin = 8;
cparamcell{1} = cparams1; cparamcell{2} = cparams2;

vert = [-1 1; 0 0];

cg = chunkgraph(vert,[1 2; 2 1],funs,cparamcell); % equivalent chunkgraph

%

% sources

ns = 10;
ts = 0.0+2*pi*rand(1,ns);
sources = [cos(ts); sin(ts)];
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 20;
ts = 0.0+2*pi*rand(1,nt);
targets = [cos(ts);sin(ts)];
targets = 0.2*targets;
targets(:,1) = [0.999999;0];
targets(:,2) = [0,0.36];
targets(:,3) = [-0.95;0];
targets(:,2) = [0,0.36];

scatter(sources(1,:),sources(2,:),'o');
scatter(targets(1,:),targets(2,:),'x');
axis equal 

chnkrtotal = merge(cg.echnks);
fkern = @(s,t) -2*chnk.helm2d.kern(zk,s,t,'D');
np = chnkrtotal.k*chnkrtotal.nch;
start = tic; D = chunkermat(chnkrtotal,fkern);
t1 = toc(start);

kerns = @(s,t) chnk.helm2d.kern(zk,s,t,'s');

% eval u on bdry

targs = chnkrtotal.r; targs = reshape(targs,2,chnkrtotal.k*chnkrtotal.nch);
targstau = tangents(chnkrtotal); 
targstau = reshape(targstau,2,chnkrtotal.k*chnkrtotal.nch);

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = targs;
kernmats = kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;


% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

fprintf('%5.2e s : time to assemble matrix\n',t1)

sysmat = D;

ndim = 1;

nsub = 40;
nch_all = horzcat(cg.echnks.nch);
npt_all = horzcat(cg.echnks.npt);
[~,nv] = size(cg.verts);
ngl = cg.k;
irowlocs = cumsum([1 cg.echnks.npt]);
rcipsav = cell(nv,1);

for ivert=1:nv
    clist = cg.vstruc{ivert}{1};
    isstart = cg.vstruc{ivert}{2};
    isstart(isstart==1) = 0;
    isstart(isstart==-1) = 1;
    nedge = length(isstart);
    iedgechunks = zeros(2,nedge);
    iedgechunks(1,:) = clist;
    iedgechunks(2,:) = 1;
    nch_use = nch_all(clist);
    iedgechunks(2,isstart==0) = nch_use(isstart==0);
    
    starind = zeros(1,2*ngl*ndim*nedge);
    corinds = cell(nedge,1);
    for i=1:nedge
        i1 = (i-1)*2*ngl*ndim+1;
        i2 = i*2*ngl*ndim;
        if(isstart(i))
        starind(i1:i2) = irowlocs(clist(i))+(1:2*ngl*ndim)-1;
            corinds{i} = 1:2*ngl;
        else
            starind(i1:i2) = irowlocs(clist(i)+1)-fliplr(0:2*ngl*ndim-1)-1;
            corinds{i} = (npt_all(clist(i))-2*ngl + 1):npt_all(clist(i));
        end
    end
    
    [Pbc,PWbc,starL,circL,starS,circS,ilist,starL1,circL1] = ...
        chnk.rcip.setup(cg.k,ndim,nedge,isstart);

    optsrcip = [];
    optsrcip.save_depth = Inf;
    tic; [R{ivert},rcipsav{ivert}] = chnk.rcip.Rcompchunk(cg.echnks,iedgechunks,fkern,ndim,cg.verts(:,ivert), ...
        Pbc,PWbc,nsub,starL,circL,starS,circS,ilist,starL1,circL1,...,
        [],[],[],[],[],optsrcip);
    toc
    
    sysmat(starind,starind) = inv(R{ivert}) - eye(2*cg.k*nedge*ndim);
end
sysmat = sysmat + eye(np);

opts = [];
opts.nsub = nsub; 
opts.rcip_savedepth = Inf;
[sysmat2,~,rcipsav] = chunkermat(cg,fkern,opts);

sysmat2 = sysmat2 + eye(np);

assert(norm(sysmat - sysmat2) < 1e-10)

[sol,flag,relres,iter] = gmres(sysmat,ubdry,np,eps*20,np);
[sol2,flag,relres,iter] = gmres(sysmat2,ubdry,np,eps*20,np);

assert(norm(sol-sol2) < 1e-10)

%

% interpolate to fine grid

ndepth = 20;
cor = cell(1,nv);
starinds = cell(1,nv);

for ivert = 1:nv
    starindtmp = rcipsav{ivert}.starind;
    starinds{ivert} = starindtmp;
    solhat = sol(starindtmp);
    [solhatinterpcell,srcinfocell,wtscell] = chnk.rcip.rhohatInterp(solhat,rcipsav{ivert},ndepth);

    rtmp = [];
    dtmp = [];
    d2tmp = [];
    ntmp = [];
    datatmp = [];
    solhatinterp = [];
    wts = [];
    for j = 1:length(srcinfocell)
        rtmp = [rtmp, srcinfocell{j}.r];
        dtmp = [dtmp, srcinfocell{j}.d];
        d2tmp = [d2tmp, srcinfocell{j}.d2];
        ntmp = [ntmp, srcinfocell{j}.n];
        if isfield(srcinfocell{j},'data')
            datatmp = [datatmp, srcinfocell{j}.data];
        end
        solhatinterp = [solhatinterp; solhatinterpcell{j}(:)];
        wts = [wts; wtscell{j}(:)];

    end
    srcinfo = struct('r',rtmp,'d',dtmp,'d2',d2tmp,'n',ntmp,'data',datatmp);

    targtemp = targets(:,:) - rcipsav{ivert}.ctr(:,1);
    targinfo = [];
    targinfo.r = targtemp;
    
    cmat = fkern(srcinfo,targinfo);
    mu = solhatinterp(:).*wts(:);
    cor{ivert} = cmat*mu;
end

%

soltemp = sol;
for ivert = 1:nv
    soltemp(starinds{ivert}) = 0;
end
    

opts = [];
opts.verb=false;
start=tic; Dsol = chunkerkerneval(chnkrtotal,fkern,soltemp,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

for ivert = 1:nv
    Dsol = Dsol + cor{ivert};
end
%

wchnkr = chnkrtotal.wts;

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkrtotal.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-10)
assert(relerr2 < 1e-10)

end


%%


%%
%----------------------------------------
%
%
% Auxiliary routines for generating boundary
%


function [r,d,d2] = circle(t)

c = cos(t(:).'); s = sin(t(:).');
r = [c; s]; 
d = [-s; c];
d2 = [-c; -s];

end
