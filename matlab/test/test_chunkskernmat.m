
%TEST_CHUNKSKERNMAT
%
% test the matrix builder

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 3;
narms = 5;
amp = 0.5;
start = tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); 
time1 = toc(start);

fprintf('%5.2e s : time to build geo\n',time1)

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

% build laplace dirichlet matrix

fkern = @(s,t,sn,tn) glapkern(s,t,sn,tn,'D');
ndims(1) = 1; ndims(2) = 1;
intparams.intorder = chunker.k;
start = tic; D = chunkskernelmat(chunker,fkern,ndims,intparams);
time1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',time1)

%

sysmat = eye(size(D))-2*D;

% built-in determinant computation

start = tic; d1 = det(sysmat); time1 = toc(start);

fprintf('%5.2e s : time for built in determinant routine\n',time1)

% using FLAM
npts = length(xs);
xhifie = chunker.chunks(1:2,:);

%
start = tic; F = rskelf(sysmat,xhifie,80,1e-10); time1 = toc(start);
fprintf('%5.2e s : time for FLAM rskelf\n',time1)
start = tic; d2 = exp(rskelf_logdet(F)); time1 = toc(start);
fprintf('%5.2e s : time for FLAM determinant\n',time1)
fprintf('\n')
fprintf('ratio of determinants : %9.6e\n',d1/d2)

