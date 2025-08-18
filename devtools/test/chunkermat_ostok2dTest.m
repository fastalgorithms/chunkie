% chunkermat_ostok2dTest0();

% function chunkermat_ostok2dTest0()

%CHUNKERMAT_OSTOK2DTEST
%
% test the matrix builder and do a basic solve

% % % % % iseed = 8675309;
% % % % % rng(iseed);
% % % % % 
% % % % % cparams = [];
% % % % % cparams.eps = 1.0e-10;
% % % % % cparams.nover = 1;
% % % % % pref = []; 
% % % % % pref.k = 16;
% % % % % narms = 3;
% % % % % amp = 0.25;
start = tic;
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];


opts = [];
opts.maxchunklen = 0.5;
pref = []; pref.nchmax = 100000;
chnkr1 = chunkerfunc(circfun, opts, pref);


chnkr2 = 0.2*chnkr1;
chnkr2 = chnkr2.reverse;

 chnkr3 = [0.5;0.0] + chnkr2;
 chnkr4 = [0.;0.5] + chnkr2;
 chnkr5 = [-0.5;0.0] + chnkr2;
 chnkr6 = [0.0;-0.5] + chnkr2;


% chnkr = merge([chnkr1, chnkr2]); % only one hole
% chnkr = merge([chnkr1, chnkr2, chnkr3]);  % only two hole
chnkr = merge([chnkr1, chnkr2, chnkr3, chnkr4, chnkr5, chnkr6]);

ns = 50;
nt = 100;

% sources = 0.05*circfun(1:ns);  % only one hole


% sources = zeros(2,2*ns); % two holes
 sources = zeros(2,5*ns);
 sources(:,1:ns) = 0.05*circfun(1:ns);
 sources(:,ns+1:2*ns) = [0.5;0.0] + 0.05*circfun(1:ns);  % stop here for two holes
 sources(:,2*ns+1:3*ns) = [-0.5;0.0] + 0.05*circfun(1:ns);
 sources(:,3*ns+1:4*ns) = [0.0;0.5] + 0.05*circfun(1:ns);
 sources(:,4*ns+1:5*ns) = [0.0;-0.5] + 0.05*circfun(1:ns);


targets = 0.90*circfun(1:nt);
% targets = targets.*repmat(rand(1,nt),2,1)*0.8;



% strengths = randn(4*ns,1); % only two holes
% strengths = randn(2*ns,1); % only one hole
 strengths = randn(10*ns,1);


% sources_n = rand(2,ns);  % only one hole
% sources_n = rand(2,2*ns); % only two holes 
sources_n = rand(2,5*ns); % five hole

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx');
plot(sources(1,:), sources(2,:), 'bo');hold on;

t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% section %% 
zk = 1.3;
kernsp = kernel('ostok', 'sp', zk);
kerns = kernel('ostok', 's', zk);
kernsink = kernel('ostok', 'sinkv', zk);

% eval u on bdry


srcinfo = []; 
srcinfo.r = sources; 
srcinfo.n = sources_n;
kernmats = kernsink.eval(srcinfo, chnkr);
% kernmats = kernsp.eval(srcinfo, chnkr);
ubdry = kernmats*strengths; 



% section %%
ubdry2 = kernsink.eval(srcinfo, chnkr)*strengths;
n1 = chnkr1.npt;
ux1 = ubdry2(1:2:2*n1);
uy1 = ubdry2(2:2:2*n1);
nx1 = chnkr1.n(1,:).';
ny1 = chnkr1.n(2,:).';

r = sum((ux1.*nx1 + uy1.*ny1).*chnkr1.wts(:))



ubdry3 = kernsink.eval(srcinfo, chnkr2)*strengths;
n2 = chnkr2.npt;
ux2 = ubdry3(1:2:2*n2);
uy2 = ubdry3(2:2:2*n2);
nx2 = chnkr2.n(1,:).';
ny2 = chnkr2.n(2,:).';

r2 = sum((ux2.*nx2 + uy2.*ny2).*chnkr2.wts(:))


ubdry4 = kernsink.eval(srcinfo, chnkr3)*strengths;
n3 = chnkr3.npt;
ux3 = ubdry4(1:2:2*n3);
uy3 = ubdry4(2:2:2*n3);
nx3 = chnkr3.n(1,:).';
ny3 = chnkr3.n(2,:).';

r3 = sum((ux3.*nx3 + uy3.*ny3).*chnkr3.wts(:))

ubdry5 = kernsink.eval(srcinfo, chnkr4)*strengths;
n4 = chnkr4.npt;
ux4 = ubdry5(1:2:2*n4);
uy4 = ubdry5(2:2:2*n4);
nx4 = chnkr4.n(1,:).';
ny4 = chnkr4.n(2,:).';

r4 = sum((ux4.*nx4 + uy4.*ny4).*chnkr4.wts(:))


ubdry6 = kernsink.eval(srcinfo, chnkr5)*strengths;
n5 = chnkr5.npt;
ux5 = ubdry6(1:2:2*n5);
uy5 = ubdry6(2:2:2*n5);
nx5 = chnkr5.n(1,:).';
ny5 = chnkr5.n(2,:).';

r5 = sum((ux5.*nx5 + uy5.*ny5).*chnkr5.wts(:))


ubdry7 = kernsink.eval(srcinfo, chnkr6)*strengths;
n6 = chnkr6.npt;
ux6 = ubdry7(1:2:2*n6);
uy6 = ubdry7(2:2:2*n6);
nx6 = chnkr6.n(1,:).';
ny6 = chnkr6.n(2,:).';

r6 = sum((ux6.*nx6 + uy6.*ny6).*chnkr6.wts(:))

% section %%
% eval u at targets

targinfo = []; 
targinfo.r = targets;
targets_n = rand(2, nt); 
targets_n = targets_n./sqrt(targets_n(1,:).^2+targets_n(2,:).^2);
targinfo.n = targets_n;
kernmatstarg = kerns.eval(srcinfo, targinfo);   %   for sprime
utarg = kernmatstarg*strengths; 
% kernmatstarg = kernsink.eval(srcinfo, targinfo);
% utarg = kernmatstarg*strengths;
%%%%% exact solution uisng analyticity not code
% section %%
% solve

% build Oscillatory Stokes dirichlet matrix
fkern = kernel('ostok', 'sp', zk);
uu = fkern.eval(chnkr, targinfo);
% section %%
start = tic; 
D = chunkermat(chnkr, fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(size(D,1)) + D;
sys = sys + normonesmat(chnkr)/sum(chnkr.wts(:));

rhs = ubdry; 
rhs = rhs(:);

start = tic; 
sol = gmres(sys,rhs,[],1e-12,1000); 
t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)


% double layer test 

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
% start=tic; 
% Dsol = chunkerkerneval(chnkr, fkern, sol, targets, opts); 
% t1 = toc(start);
% fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)
% 
% relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
% fprintf('relative frobenius error %5.2e\n',relerr);
% 
% assert(relerr < 1e-10);

%%%%%%%% DOUBLE LAYER TEST

% fkernd =  kernel('ostok', 'd', zk);
% Dsol = chunkerkerneval(chnkr, fkernd, sol, targets, opts);
% relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
% fprintf('relative frobenius error %5.2e\n',relerr);
% 
% fprintf('relative frobenius error %5.2e\n',relerr);

%%%%%%%   SPRIME TEST

fkerns =  kernel('ostok', 's', zk);
Dsol = chunkerkerneval(chnkr, fkerns, sol, targets, opts);
relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('relative frobenius error %5.2e\n',relerr);

fprintf('relative frobenius error %5.2e\n',relerr);

% assert(relerr < 1e-10);

% end