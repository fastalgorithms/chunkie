%CHUNKERKERNEVAL_MULTIDENSTEST test chunkerkerneval with several densities
% at once, against one density at a time

chunkerkerneval_multidensTest0();


function chunkerkerneval_multidensTest0()

rng(1);
cparams = [];
cparams.eps = 1e-6;
pref = [];
pref.k = 16;
chnkr = chunkerfunc(@(t) starfish(t,3,0.2),cparams,pref);

% a few targets near and far. there are enough sources for the fmm to be
% used
nt = 40;
th = 2*pi*rand(1,nt);
rr = 0.5 + 2*rand(1,nt);
targs = [];
targs.r = [rr.*cos(th); rr.*sin(th)];
targs.n = [cos(th); sin(th)];
targs.d = [-sin(th); cos(th)];

ndens = 2;
zk = 2.1;
gh = @(s,t) chnk.helm2d.green(zk,s,t);
fmmh = @(eps,s,pg,t,pgt) hfmm2d(eps,zk,s,pg,t,pgt);
kerns = {kernel('lap','sgrad'), kernel('helm','c',zk,[1 1i]), ...
    kernel('helm','all',zk,[1 2; 3 4i]), kernel('stok','strac',1.3), ...
    kernel([kernel('lap','s'); kernel('helm','d',zk)]), ...
    [1 2; 3 4]*kernel('stok','svel',1.3), ...
    kernel('stok','dvel',1.3)*@(s) reshape(s.r,2,1,[]).*[1 1], ...
    kernel.green(gh,'trans_rep',[1 2i],'fmm',fmmh)};

for i = 1:numel(kerns)
    K = kerns{i};
    dens = randn(K.opdims(2)*chnkr.npt,ndens);
    u = chunkerkerneval(chnkr,K,dens,targs);
    assert(isequal(size(u),[K.opdims(1)*nt ndens]),'multidens kernel %d: size',i);
    for l = 1:ndens
        u1 = chunkerkerneval(chnkr,K,dens(:,l),targs);
        err = norm(u(:,l) - u1)/norm(u1);
        assert(err < 1e-10,'multidens kernel %d: %.2e',i,err);
    end
end

% the other quadratures, and the density as opdims(2) x k x nch x ndens
K = kernel('helm','s',zk);
dens = randn(chnkr.k,chnkr.nch,ndens);
tnear = struct('r',targs.r(:,rr < 1));
tnear.r = tnear.r(:,1:min(4,end));
optss = {struct('accel',false), struct('forcesmooth',true), ...
    struct('forceadap',true), struct('forcepquad',true,'side','i'), ...
    struct('flam',true)};
for j = 1:numel(optss)
    u = chunkerkerneval(chnkr,K,reshape(dens,1,chnkr.k,chnkr.nch,ndens), ...
        tnear,optss{j});
    for l = 1:ndens
        u1 = chunkerkerneval(chnkr,K,reshape(dens(:,:,l),[],1),tnear,optss{j});
        assert(norm(u(:,l) - u1) < 1e-10*norm(u1),'multidens opts %d',j);
    end
end
cormat = chunkerkernevalmat(chnkr,K,targs,struct('corrections',true));
u = chunkerkerneval(chnkr,K,reshape(dens,[],ndens),targs,struct('cormat',cormat));
for l = 1:ndens
    u1 = chunkerkerneval(chnkr,K,reshape(dens(:,:,l),[],1),targs, ...
        struct('cormat',cormat));
    assert(norm(u(:,l) - u1) < 1e-10*norm(u1),'multidens cormat');
end

% chunkgraph with a kernel per region
verts = [1 0 -1 0; 0 1 0 -1];
cg = chunkgraph(verts,[1:4; [2:4 1]]);
Kr = [kernel('lap','s'); kernel('helm','s',zk)];
tg = [0.1 3 -0.2 2.5; 0.05 0.3 -0.1 -2];
dens = randn(cg.npt,ndens);
u = chunkerkerneval(cg,Kr,dens,tg);
for l = 1:ndens
    u1 = chunkerkerneval(cg,Kr,dens(:,l),tg);
    assert(norm(u(:,l) - u1) < 1e-10*norm(u1),'multidens chunkgraph');
end

% the fmms against direct evaluation
src = []; th = 2*pi*rand(1,200);
src.r = [cos(th); sin(th)];  src.n = src.r;
targ = []; targ.r = 3*rand(2,150) - 1.5;  targ.n = randn(2,150);
targ.n = targ.n./vecnorm(targ.n);  targ.d = [-targ.n(2,:); targ.n(1,:)];
kerns = {};
for type = {'all','trep','trep_p','c2trans','s','d','sp','dp','c','cp'}
    switch type{1}
        case {'all','c2trans'}, a = {[1 2i; 3 4]};
        case {'trep','trep_p','c','cp'}, a = {[1 2i]};
        otherwise, a = {};
    end
    kerns{end+1} = kernel('helm',type{1},zk,a{:}); %#ok<AGROW>
end
for type = {'s','d','sp','dp','sgrad','dgrad','stau'}
    kerns{end+1} = kernel('lap',type{1}); %#ok<AGROW>
end
for type = {'svel','spres','strac','sgrad','dvel','dpres','dtrac','dgrad'}
    kerns{end+1} = kernel('stok',type{1},1.3); %#ok<AGROW>
end
kerns{end+1} = [1 2; 3 4]*kernel('stok','svel',1.3);
kerns{end+1} = kernel('stok','dvel',1.3)*[2 1; 0 1];
for i = 1:numel(kerns)
    K = kerns{i};
    A = K.eval(src,targ);
    sigma = randn(size(A,2),ndens);
    u = K.fmm(1e-12,src,targ,sigma);
    assert(norm(u - A*sigma,'fro') < 1e-9*norm(A*sigma,'fro'), ...
        'fmm kernel %d',i);
end

% the gradient and hessian outputs of the helmholtz fmm
sigma = randn(size(src.r,2),1);
[pot,grad,hess] = chnk.helm2d.fmm(1e-12,zk,src,targ,'s',sigma);
[v,g,h] = chnk.helm2d.green(zk,src.r,targ.r);
ref = {v*sigma, [g(:,:,1)*sigma g(:,:,2)*sigma].', ...
    [h(:,:,1)*sigma h(:,:,2)*sigma h(:,:,3)*sigma].'};
out = {pot,grad,hess};
for j = 1:3
    assert(norm(out{j}(:) - ref{j}(:)) < 1e-9*norm(ref{j}(:)), ...
        'helm fmm output %d',j);
end

% the pressure and gradient outputs of the stokes fmm, several densities
mu = 1.3;
sigma = randn(2*size(src.r,2),ndens);
[vel,pres,grad] = chnk.stok2d.fmm(1e-12,mu,src,targ,'svel',sigma);
ref = {kernel('stok','svel',mu).eval(src,targ)*sigma, ...
    kernel('stok','spres',mu).eval(src,targ)*sigma, ...
    kernel('stok','sgrad',mu).eval(src,targ)*sigma};
out = {vel,pres,reshape(grad,[],ndens)};
for j = 1:3
    assert(norm(out{j} - ref{j},'fro') < 1e-9*norm(ref{j},'fro'), ...
        'stok fmm output %d',j);
end

end
