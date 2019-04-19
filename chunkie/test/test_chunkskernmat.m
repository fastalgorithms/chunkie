
%TEST_CHUNKSKERNMAT
%
% test the matrix builder

addpaths_loc();

cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 10;
amp = 0.5;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

figure(1)
plot(chnkr)
hold on
quiver(chnkr)
axis equal

fprintf('%5.2e s : time to build geo\n',t1)

% build laplace dirichlet matrix

fkern = @(s,t,stau,ttau) glapkern(s,t,stau,ttau,'D');
opdims(1) = 1; opdims(2) = 1;
intparams.intorder = chnkr.k;
start = tic; D = chunkskernmat(chnkr,fkern,opdims,intparams);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + D;

rhs = cos(sum(abs(chnkr.r).^2,1)); rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-12,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)