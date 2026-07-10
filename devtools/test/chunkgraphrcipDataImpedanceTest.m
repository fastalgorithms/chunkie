chunkgraphrcip_Data_Impedance_Test();


function chunkgraphrcip_Data_Impedance_Test()

%CHUNKGRAPHRCIP_DATA_IMPEDANCE_TEST
%
% Solving the interior Helmholtz problem with impedance BCs:
%
%       alpha(x) du/dn + i beta(x) u = 0 , 
%
% where alpha and beta are allowed to vary along the boundary. We solve
% this problem using the following scaled Laplace single layer ansatz: 
%   
%           u = 2 \int G(x,y) \mu(y) / alpha(y) ds(y) , 
% 
% where G(x,y) is the Helmholtz Green's function and alpha is bounded 
% below by a positive number.
%
% This problem is solved on a square chunkgraph where alpha and beta are 
% stored as chunkgraph data. This test checks whether chunkgraph data
% is properly added and stored and whether RCIP is properly interpolating 
% source and target data.

rng(1);

t0 = tic;
vertsq = [1 1 -1 -1; -1 1 1 -1];
edges = [1 2 3 4; 2 3 4 1];

cparams = [];
cparams.eps = 1.0e-10;
[cgrph] = chunkgraph(vertsq,edges,[],cparams);

cgrph = makedatarows(cgrph,2);

% Impedance coefficients, chosen arbitrarily
cgrph.data(1,:) = [cgrph.echnks(1).r(2,:).^2+2, ...
    cgrph.echnks(2).r(1,:).^2+2, ...
    cgrph.echnks(3).r(2,:).^2+2, ...
    cgrph.echnks(4).r(1,:).^2+2];  % alpha 

cgrph.data(2,:) = [cgrph.echnks(1).r(1,:), ...
    cgrph.echnks(2).r(2,:)+2, ...
    cgrph.echnks(3).r(2,:)-2, ...
    cgrph.echnks(4).r(1,:).^2]; % beta

t1 = toc(t0);

fprintf('%5.2e s : time to build geo\n',t1)

fprintf('Number of points: %d\n',cgrph.npt)
fprintf('Number of edges: %d\n',cgrph.nedge)

zk = 1.1;

% sources

sources = [[-2.3;0],[1.6;1],[3.1;-0.7]];
strengths = rand(size(sources,2),1);

% targets

nt = 3;
ts = 0.0+2*pi*rand(nt,1);
targets = 2*rand(2,nt)-1;

% plot geo and sources

xs = cgrph.r(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = cgrph.r(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));

figure(1)
clf
hold off
plot(cgrph)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

kern = @(s,t) kerntot(s,t,zk);
kerns = @(s,t) 2*chnk.helm2d.kern(zk,s,t,'s');
kernsp = @(s,t) 2*chnk.helm2d.kern(zk,s,t,'sp');
kernrep = @(s,t) 2*chnk.helm2d.kern(zk,s,t,'s') ./ s.data(1,:);

% eval u on bdry

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = cgrph.r; targinfo.n = cgrph.n; targinfo.data = cgrph.data;
kernmats = cgrph.data(1,:).'.*kernsp(srcinfo,targinfo) ...
    + 1i*cgrph.data(2,:).'.*kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

%

% build laplace dirichlet matrix

opts = []; opts.rcip = true; opts.sing = 'log'; opts.nsub_or_tol = 40;
start = tic; 
A = chunkermat(cgrph,kern,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = eye(cgrph.npt) + A;

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare

opts = [];
start=tic; Dsol = chunkerkerneval(cgrph,kernrep,sol2,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wcgrph = cgrph.wts;

relerr = norm(utarg-Dsol,'fro')/norm(utarg,'fro');
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wcgrph(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-10);

end


function val = kerntot(s,t,zk)

alphas = s.data(1,:);
alphat = t.data(1,:);

beta = t.data(2,:);

val = 2*alphat(:).*chnk.helm2d.kern(zk,s,t,'sp') ./ alphas ...
    + 2i*beta(:).*chnk.helm2d.kern(zk,s,t,'s') ./ alphas;

end