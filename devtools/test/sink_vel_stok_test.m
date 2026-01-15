start = tic;
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];


opts = [];
opts.maxchunklen = 0.5;
pref = []; pref.nchmax = 100000;
chnkr1 = chunkerfunc(circfun, opts, pref);

% % % % chnkr2 = 0.05*chnkr1;
% % % % chnkr2 = chnkr2.reverse;

chnkr = chnkr1;
% % % % chnkr = merge([chnkr1, chnkr2]);

ns = 10;
nt = 100;

sources = 1.3*circfun(1:ns);

targets = 0.75*circfun(1:nt);


strengths = randn(ns,1);
sources_n = rand(2,ns);
plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx');hold on;
plot(sources(1,:), sources(2,:), 'bo');hold on;

t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)
 % % % % % 
 % % % % % plot(chnkr, 'r.'); hold on;
 % % % % % plot(targets(1,:), targets(2,:), 'kx');
 % % % % % plot(sources(1,:), sources(2,:), 'bo');
 % % % % % axis equal

% section %% 
zk = 1.3;

kernsink = kernel('ostok', 'sinkv', zk);

% eval u on bdry


srcinfo = []; 
srcinfo.r = sources; 
srcinfo.n = sources_n;
kernmats = kernsink.eval(srcinfo, chnkr);
ubdry = kernmats*strengths; 


% eval u at targets

kernmatstarg = kernsink.eval(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

% build Stokes dirichlet matrix

start = tic; 
D = chunkermat(chnkr,kernsink);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(size(D,1)) + D;

sys = sys + normonesmat(chnkr)/sum(chnkr.wts(:));

rhs = ubdry; 
rhs = rhs(:);
start = tic; 
sol = gmres(sys,rhs,[],1e-14,100); 
t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)
