quasiperiodicTest0();
quasiperiodicTest1();
quasiperiodicTest2();
quasiperiodicTest_flex_clamped();
quasiperiodicTest_flex_free();
quasiperiodicTest_flex_supported();

function quasiperiodicTest0()
% test the representation

% problem parameters
d= 1;
zk = 1;
kappa = [pi-0.1-1i,0.2,-0.1,0.25i];
coefs = [1;0.5];
coefa = [1,0.3;-1i,0.2];
nkappa = length(kappa);

src = []; src.r = [[0;-1.1],[1;-1],[0.1;-0.3]]; src.n = [[1;-2],[2;-1],[1;0]];
targ = []; targ.r = [[1.1;0.3],[2;0]]; targ.n = [[-1;0.3],[0.1;1]];
ns = size(src.r,2);
nt = size(targ.r,2);


% kernels
skern = kernel('hq','s',zk,kappa,d);
dkern = kernel('hq','d',zk,kappa,d);
ckern = kernel('hq','c',zk,kappa,d,coefs);

spkern = kernel('hq','sp',zk,kappa,d);
dpkern = kernel('hq','dp',zk,kappa,d);
cpkern = kernel('hq','cp',zk,kappa,d,coefs);

quas_param = skern.params.quas_param;

allkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'all',quas_param,coefa));
trkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'trans_rep',quas_param,coefs));
trpkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'trans_rep_prime',quas_param,coefs));

sgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'sgrad',quas_param));
dgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'dgrad',quas_param));
cgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'cgrad',quas_param,coefs));
trgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'trans_rep_grad',quas_param,coefs));

c2trkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'c2trans',quas_param,coefs));

% evaluate kernels
sval = skern.eval(src,targ);
dval = dkern.eval(src,targ);
cval = ckern.eval(src,targ);
trval = trkern.eval(src,targ);

spval = spkern.eval(src,targ);
dpval = dpkern.eval(src,targ);
cpval = cpkern.eval(src,targ);
trpval = trpkern.eval(src,targ);

sgval = sgkern.eval(src,targ);
dgval = dgkern.eval(src,targ);
cgval = cgkern.eval(src,targ);
trgval = trgkern.eval(src,targ);

allval = allkern.eval(src,targ);

c2trval = c2trkern.eval(src,targ);


% test quasiperiodicitiy
targ_s = targ; targ_s.r = targ.r + [d;0];
svalshift = skern.eval(src,targ_s);
assert(norm(svalshift - repmat(exp(1i*kappa(:)*d),nt,1).*sval) <  1e-12)

% test combined
assert(norm((coefs(1)*dval + coefs(2)*sval)-cval) <  1e-12)
assert(norm((coefs(1)*dpval + coefs(2)*spval)-cpval) <  1e-12)

% test all
all_assemb = zeros(2,nkappa*nt,2*ns);

all_assemb(1,:, 1:2:end) = coefa(1,1)*dval;
all_assemb(1,:, 2:2:end) = coefa(1,2)*sval;
all_assemb(2,:, 1:2:end) = coefa(2,1)*dpval;
all_assemb(2,:, 2:2:end) = coefa(2,2)*spval;
all_assemb = reshape(all_assemb,[2,nkappa,nt,2*ns]);
all_assemb = permute(all_assemb,[2,1,3,4]);
all_assemb = reshape(all_assemb, size(allval) );

assert(norm(all_assemb - allval) <  1e-12)

% test transmission representation
tr_assemb = zeros(nkappa*nt,2,ns);
tr_assemb(:,1,:) = coefs(1)*dval;
tr_assemb(:,2,:) = coefs(2)*sval;
assert(norm(tr_assemb(:,:)-trval) <  1e-12)
trp_assemb = zeros(nkappa*nt,2,ns);
trp_assemb(:,1,:) = coefs(1)*dpval;
trp_assemb(:,2,:) = coefs(2)*spval;
assert(norm(trp_assemb(:,:)-trpval) <  1e-12)

% test gradients
assert(norm((coefs(1)*dgval + coefs(2)*sgval)-cgval) <  1e-12)
trg_assemb = zeros(2*nkappa*nt,2,ns);
trg_assemb(:,1,:) = coefs(1)*dgval;
trg_assemb(:,2,:) = coefs(2)*sgval;
assert(norm(trg_assemb(:,:)-trgval) <  1e-12)

% test c2trans
c2trval_assemb = zeros(nkappa,2,nt,ns);
c2trval_assemb(:,1,:,:) = reshape(coefs(1)*dval + coefs(2)*sval,nkappa,1,nt,ns);
c2trval_assemb(:,2,:,:) = reshape(coefs(1)*dpval + coefs(2)*spval,nkappa,1,nt,ns);
c2trval_assemb = reshape(c2trval_assemb,[],ns);
assert(norm(c2trval_assemb-c2trval) <  1e-12)

% test ising parameter
skern1a = kernel('hq','s',zk,kappa,d,[],[],0);
skern1b = kernel('h','s',zk);
sval1b = skern1b.eval(src,targ); 
sval1b = repmat(reshape(sval1b,1,nt,ns), nkappa,1,1);
sval1 = skern1a.eval(src,targ) + reshape(sval1b,nt*nkappa,ns);
assert(norm(sval-sval1) <  1e-12)

dkern1a = kernel('hq','d',zk,kappa,d,[],[],0);
dkern1b = kernel('h','d',zk);
dval1b = dkern1b.eval(src,targ); 
dval1b = repmat(reshape(dval1b,1,nt,ns), nkappa,1,1);
dval1 = dkern1a.eval(src,targ) + reshape(dval1b,nt*nkappa,ns);
assert(norm(dval-dval1) <  1e-12)

allkern1a = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'all',quas_param,coefa,0));
allkern1b = kernel(@(s,t) chnk.helm2d.kern(zk,s,t,'all',coefa));
allval1a = allkern1a.eval(src,targ);
allval1b = allkern1b.eval(src,targ); 
allval1b = reshape(allval1b,1,2*nt,2*ns); allval1b = repmat(allval1b,length(kappa),1,1);
allval1b = reshape(allval1b,length(kappa)*2*nt,2*ns);
assert(norm(allval1a+allval1b -allval) <  1e-12)
end


function quasiperiodicTest1()
% test an integral equation


% problem parameters
d= 8;
zk = 1;
kappa = .05+.1i;

% setup geometry
nch = 2^3;
A = .2;
cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
chnkr = chunkerfuncuni(@(t) sin_func(t,d,A),nch,cparams);
chnkr = reverse(chnkr);

% setup system
skern = kernel('hq','s',zk,kappa,d);
dkern = kernel('hq','d',zk,kappa,d);

sysmat = chunkermat(chnkr,dkern);
sysmat = sysmat + .5*eye(size(sysmat,2));

% solve
src = struct("r",[3*d/4;-1.5]);

rhs = skern.eval(src,chnkr);

soln = sysmat\rhs;

% check analytic solution
targ = [10*d/3;2];
utrue = skern.eval(src,struct("r",targ));
u = chunkerkerneval(chnkr,dkern,soln,targ);

assert(abs(u-utrue) < 1e-10)

% check gradient
quas_param = skern.params.quas_param;
sgkern = @(s,t) chnk.helm2dquas.kern(zk,s,t,'sgrad',quas_param);
dgkern = @(s,t) chnk.helm2dquas.kern(zk,s,t,'dgrad',quas_param);

ugradtrue = sgkern(src,struct("r",targ));
ugrad = chunkerkerneval(chnkr,dgkern,soln,targ);

assert(norm(ugrad-ugradtrue) < 1e-8)

opts = [];
opts.forcepquad = true; opts.side = 'e';
t = 0;
targp = sin_func(t,d,A)+1e-2;
uscatp = chunkerkerneval(chnkr,dkern,soln,targp,opts);
uinp = skern.eval(src,struct("r",targp));
assert(abs(uinp-uscatp) < 1e-10)

% ploting
Lplot = d/1.5;
nplot = 50;
xts = linspace(-Lplot,Lplot,nplot); yts=xts+Lplot/2;
[X,Y] = meshgrid(xts,yts);
targ = [X(:).';Y(:).'];

% targ = [2+0*yts;yts];
ntarg = size(targ,2);
tic;
yc = sin_func(X(:),d,A); yc = yc(2,:);
iup = Y(:).'>yc;
toc
targup = targ(:,iup.');

tic;
uin = NaN*zeros(ntarg,1)+NaN*1i;
src.n = [0;1];
uin(iup) = skern.eval(src,struct("r",targup));
uscat = NaN*zeros(ntarg,1)+NaN*1i;
opts = []; opts.forcesmooth = false; opts.eps = 1e-6;
opts.forcepquad = true; opts.side = 'e';
uscat(iup) = chunkerkerneval(chnkr,dkern,soln,targup,opts);
toc
%
utot = uin-uscat;
maxu = max(abs(uin));
figure(1);clf
subplot(1,3,1)
h = pcolor(X,Y,reshape(imag(uin),nplot,[])); set(h,'edgecolor','none')
title('$u_{\rm{true}}$','Interpreter','latex'); set(gca,'fontsize',14)
colorbar; clim([-maxu,maxu])
hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off

subplot(1,3,2)
h= pcolor(X,Y,reshape(imag(uscat),nplot,[])); set(h,'edgecolor','none')
title('$u$','Interpreter','latex'); set(gca,'fontsize',14)
colorbar; clim([-maxu,maxu])
hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off

subplot(1,3,3)
h = pcolor(X,Y,reshape(log10(abs(utot)),nplot,[])); set(h,'edgecolor','none')
title('$\log_{10} $ error','Interpreter','latex'); set(gca,'fontsize',14)
colorbar; 
hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off


end

function quasiperiodicTest2()
% test an integral equation for the quasiperiodic Laplace problem

% problem parameters
d = 8;
kappa = .05-.1i;

% setup geometry
nch = 2^3;
A = .2;
cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
chnkr = chunkerfuncuni(@(t) sin_func(t,d,A),nch,cparams);
chnkr = reverse(chnkr);

% chnkr = chunkerfunc(@starfish);

% setup system
skern = kernel('lq','s',kappa,d);
spkern = kernel('lq','sp',kappa,d);

sysmat = chunkermat(chnkr,spkern);
sysmat = sysmat - .5*eye(size(sysmat,2));

% solve
src = struct("r",[3*d/4;-1.5]);

rhs = spkern.eval(src,chnkr);

soln = sysmat\rhs;

% check analytic solution
targ = [10*d/3;2];
utrue = skern.eval(src,struct("r",targ));
u = chunkerkerneval(chnkr,skern,soln,targ);

assert(abs(u-utrue) < 1e-10)

% % ploting
% Lplot = d/1.5;
% nplot = 50;
% xts = linspace(-Lplot,Lplot,nplot); yts=xts+Lplot/2;
% [X,Y] = meshgrid(xts,yts);
% targ = [X(:).';Y(:).'];
% 
% % targ = [2+0*yts;yts];
% ntarg = size(targ,2);
% tic;
% yc = sin_func(X(:),d,A); yc = yc(2,:);
% iup = Y(:).'>yc;
% % iup = ~chunkerinterior(chnkr,targ);
% toc
% targup = targ(:,iup.');
% 
% tic;
% uin = NaN*zeros(ntarg,1)+NaN*1i;
% src.n = [0;1];
% uin(iup) = skern.eval(src,struct("r",targup));
% uscat = NaN*zeros(ntarg,1)+NaN*1i;
% opts = []; opts.forcesmooth = false; opts.eps = 1e-6;
% uscat(iup) = chunkerkerneval(chnkr,skern,soln,targup,opts);
% toc
% %
% utot = uin-uscat;
% maxu = max(abs(uin));
% figure(2);clf
% subplot(1,3,1)
% h = pcolor(X,Y,reshape(imag(uin),nplot,[])); set(h,'edgecolor','none')
% title('$u_{\rm{true}}$','Interpreter','latex'); set(gca,'fontsize',14)
% colorbar; clim([-maxu,maxu])
% hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off
% 
% subplot(1,3,2)
% h= pcolor(X,Y,reshape(imag(uscat),nplot,[])); set(h,'edgecolor','none')
% title('$u$','Interpreter','latex'); set(gca,'fontsize',14)
% colorbar; clim([-maxu,maxu])
% hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off
% 
% subplot(1,3,3)
% h = pcolor(X,Y,reshape(log10(abs(utot)),nplot,[])); set(h,'edgecolor','none')
% title('$\log_{10} $ error','Interpreter','latex'); set(gca,'fontsize',14)
% colorbar; 
% hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off


end

function quasiperiodicTest_flex_clamped()
% test an integral equation for the quasiperiodic clamped plate problem

% problem parameters
d = 2; zk = 7; kappa = 0.2-0.1*1i;
l = 2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).', zk, d, kappa, exp(1i*kappa*d), a, M, l+1);

% setup geometry
nch = 2^3; A = .2;
cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
chnkr = reverse(chunkerfuncuni(@(t) sin_func(t,d,A), nch, cparams));
curv = signed_curvature(chnkr); curv = curv(:);

opts_s = []; opts_s.sing = 'smooth'; opts_s.quad = 'native';
opts_l = []; opts_l.sing = 'log';

% build system: quasi periodic part uses smooth quadrature (ising=0) since
% the self-interaction and nsub nearest periodic copies are subtracted and
% handled separately via the free-space kernel with log quadrature
ising = 0; nsub = 1;
alpha = exp(1i*kappa*d);

fkern   = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate', kappa, d, sn, [], [], l, ising, nsub);
fkern_0 = @(s,t) chnk.flex2d.kern(    zk, s, t, 'clamped_plate');

sys   = reshape(chunkermat(chnkr, fkern,   opts_s), 1, 2*chnkr.npt, 2*chnkr.npt);
sys_0 = reshape(chunkermat(chnkr, fkern_0, opts_l), 1, 2*chnkr.npt, 2*chnkr.npt);
for ii = [-nsub:-1, 1:nsub]
    sys_0 = sys_0 + alpha^ii .* reshape(chunkerkernevalmat(chnkr + ii*[d;0], fkern_0, chnkr), 1, 2*chnkr.npt, 2*chnkr.npt);
end

sys = sys + sys_0 - reshape(0.5*eye(2*chnkr.npt), 1, 2*chnkr.npt, 2*chnkr.npt);
sys(:,2:2:end,1:2:end) = sys(:,2:2:end,1:2:end) + reshape(curv .* eye(chnkr.npt), 1, chnkr.npt, chnkr.npt);

% solve
ising = 1;
src    = struct('r', [0; -1.5]);
bskern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_bcs', kappa, d, sn, [], [], l, ising);
skern  = @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',                  kappa, d, sn, [], [], l, ising);

rhs = -bskern(src, chnkr);
sol = squeeze(sys) \ rhs;

% check: total field vanishes exterior => uscat + uin = 0
targ    = struct('r', [3*d/4; 1.5]);
ikern   = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval', kappa, d, sn, [], [], l, 0);
ikern_0 = @(s,t) chnk.flex2d.kern(    zk, s, t, 'clamped_plate_eval');

wts    = repmat(chnkr.wts(:).', 2, 1);
uscat  = (ikern(chnkr, targ) .* wts(:).' + chunkerkernevalmat(chnkr, ikern_0, targ, opts_l)) * sol;
uin    = skern(src, targ);
assert(abs(uscat + uin) < 1e-8)

end


function quasiperiodicTest_flex_free()
% test an integral equation for the quasiperiodic free plate problem

% problem parameters
d = 2; zk = 7; kappa = 0.2-0.1*1i; nu = 0.3;
l = 2; N = 40; a = 15; M = 1e4;
sn            = chnk.flex2dquas.latticecoefs((0:N).', zk, d, kappa, exp(1i*kappa*d), a, M, l+1);
[s0_l, sn_l] = chnk.lap2dquas.latticecoefs((1:N), d, kappa, l);

% setup geometry
nch = 2^3; A = .2;
cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
chnkr = reverse(chunkerfuncuni(@(t) sin_func(t,d,A), nch, cparams));

opts_s = []; opts_s.sing = 'smooth'; opts_s.quad = 'native';
opts_l = []; opts_l.sing = 'log';
opts_pv = []; opts_pv.sing = 'pv';

% build system: quasi periodic part uses smooth quadrature (ising=0) since
% the self-interaction and nsub nearest periodic copies are subtracted and
% handled separately via the free-space kernels with log/pv quadrature
ising = 0; nsub = 1;
alpha = exp(1i*kappa*d);

fkern   = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate', kappa, d, sn, s0_l, sn_l, l, ising, nu, nsub);
double  = @(s,t) chnk.lap2dquas.kern(      s, t, 'd',          kappa, d, s0_l, sn_l, l, ising, nsub);
hilbert = @(s,t) chnk.lap2dquas.kern(      s, t, 'hilb',       kappa, d, s0_l, sn_l, l, ising, nsub);
fkern_0   = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);
double_0  = @(s,t) chnk.lap2d.kern(      s, t, 'd');
hilbert_0 = @(s,t) chnk.lap2d.kern(      s, t, 'hilb');

sysmat1 = reshape(chunkermat(chnkr, fkern,   opts_s), 1, 4*chnkr.npt, 2*chnkr.npt);
D       = reshape(chunkermat(chnkr, double,  opts_s), 1, chnkr.npt, chnkr.npt);
H       = reshape(chunkermat(chnkr, hilbert, opts_s), 1, chnkr.npt, chnkr.npt);

sysmat1_0 = reshape(chunkermat(chnkr, fkern_0,   opts_l),  1, 4*chnkr.npt, 2*chnkr.npt);
D_0       = reshape(chunkermat(chnkr, double_0,  opts_l),  1, chnkr.npt, chnkr.npt);
H_0       = reshape(chunkermat(chnkr, hilbert_0, opts_pv), 1, chnkr.npt, chnkr.npt);
for ii = [-nsub:-1, 1:nsub]
    sysmat1_0 = sysmat1_0 + alpha^ii .* reshape(chunkerkernevalmat(chnkr + ii*[d;0], fkern_0,   chnkr), 1, 4*chnkr.npt, 2*chnkr.npt);
    D_0       = D_0       + alpha^ii .* reshape(chunkerkernevalmat(chnkr + ii*[d;0], double_0,  chnkr), 1, chnkr.npt, chnkr.npt);
    H_0       = H_0       + alpha^ii .* reshape(chunkerkernevalmat(chnkr + ii*[d;0], hilbert_0, chnkr), 1, chnkr.npt, chnkr.npt);
end

sysmat1 = sysmat1 + sysmat1_0;
D = D + D_0;
H = H + H_0;

s11b   = permute(sysmat1(:,3:4:end,1:2:end), [2,3,1]);
s21b   = permute(sysmat1(:,4:4:end,1:2:end), [2,3,1]);
D_perm = permute(D, [2,3,1]);
H_perm = permute(H, [2,3,1]);

k11tmp = permute(pagemtimes(s11b, H_perm) - 2*((1+nu)/2)^2*pagemtimes(D_perm, D_perm), [3,1,2]);
k21tmp = permute(pagemtimes(s21b, H_perm), [3,1,2]);

sysmat = zeros(1, 2*chnkr.npt, 2*chnkr.npt);
sysmat(:,1:2:end,1:2:end) = sysmat1(:,1:4:end,1:2:end) + k11tmp;
sysmat(:,2:2:end,1:2:end) = sysmat1(:,2:4:end,1:2:end) + k21tmp;
sysmat(:,1:2:end,2:2:end) = sysmat1(:,1:4:end,2:2:end) + sysmat1(:,3:4:end,2:2:end);
sysmat(:,2:2:end,2:2:end) = sysmat1(:,2:4:end,2:2:end) + sysmat1(:,4:4:end,2:2:end);

Djump = reshape(kron(eye(chnkr.npt), [-1/2 + (1/8)*(1+nu)^2, 0; 0, 1/2]), 1, 2*chnkr.npt, 2*chnkr.npt);
sys = Djump + sysmat;

% solve
ising = 1;
src    = struct('r', [0; -1.5]);
bskern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs', kappa, d, sn, s0_l, sn_l, l, ising, nu);
skern  = @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',              kappa, d, sn, s0_l, sn_l, l, ising);

rhs = -bskern(src, chnkr);
sol = squeeze(sys) \ rhs;

% build combined density [phi; H*phi; psi] for free_plate_eval
dens_comb = zeros(3*chnkr.npt, size(rhs,2));
dens_comb(1:3:end,:) = sol(1:2:end,:);
dens_comb(2:3:end,:) = squeeze(H) * sol(1:2:end,:);
dens_comb(3:3:end,:) = sol(2:2:end,:);

% check: total field vanishes exterior => uscat + uin = 0
targ    = struct('r', [3*d/4; 1.5]);
ikern   = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval', kappa, d, sn, s0_l, sn_l, l, 0, nu);
ikern_0 = @(s,t) chnk.flex2d.kern(    zk, s, t, 'free_plate_eval', nu);

wts   = repmat(chnkr.wts(:).', 3, 1);
uscat = (ikern(chnkr, targ) .* wts(:).' + chunkerkernevalmat(chnkr, ikern_0, targ, opts_l)) * dens_comb;
uin   = skern(src, targ);
assert(abs(uscat + uin) < 1e-8)

end


function quasiperiodicTest_flex_supported()
% test an integral equation for the quasiperiodic supported plate problem

% problem parameters
d = 2; zk = 7; kappa = 0.2-0.1*1i; nu = 0.3;
l = 2; N = 40; a = 15; M = 1e4;
sn            = chnk.flex2dquas.latticecoefs((0:N).', zk, d, kappa, exp(1i*kappa*d), a, M, l+1);
[s0_l, sn_l] = chnk.lap2dquas.latticecoefs((1:N), d, kappa, l);

% setup geometry (curvature arc-length derivatives stored in chnkr.data)
nch = 2^3; A = .2;
cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
chnkr = reverse(chunkerfuncuni(@(t) sin_func(t,d,A), nch, cparams));
chnkr = makedatarows(chnkr, 2);
curv  = signed_curvature(chnkr); curv = curv(:);
chnkr.data(1,:,:) = arclengthder(chnkr, curv);
chnkr.data(2,:,:) = arclengthder(chnkr, squeeze(chnkr.data(1,:,:)));

opts_s = []; opts_s.sing = 'smooth'; opts_s.quad = 'native';
opts_l = []; opts_l.sing = 'log';

% build system: quasi periodic part uses smooth quadrature (ising=0) since
% the self-interaction and nsub nearest periodic copies are subtracted and
% handled separately via the free-space kernels with log/smooth quadrature
ising = 0; nsub = 1;
alpha = exp(1i*kappa*d);
c0    = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

fkern    = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate',        kappa, d, sn, [], [], l, ising, nu, nsub);
fkern_0l = @(s,t) chnk.flex2d.kern(    zk, s, t, 'supported_plate_log',    nu);
fkern_0s = @(s,t) chnk.flex2d.kern(    zk, s, t, 'supported_plate_smooth', nu);

M3 = reshape(chunkermat(chnkr, fkern, opts_s), 1, 2*chnkr.npt, 2*chnkr.npt);

M  = chunkermat(chnkr, fkern_0l, opts_l);
M2 = reshape(chunkermat(chnkr, fkern_0s, opts_s), 1, chnkr.npt, chnkr.npt);
M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + c0 .* curv.^2 .* eye(chnkr.npt);
M  = reshape(M - 0.5*eye(2*chnkr.npt), 1, 2*chnkr.npt, 2*chnkr.npt);
for ii = [-nsub:-1, 1:nsub]
    M  = M  + alpha^ii .* reshape(chunkerkernevalmat(chnkr + ii*[d;0], fkern_0l, chnkr),                            1, 2*chnkr.npt, 2*chnkr.npt);
    M2 = M2 + alpha^ii .* reshape(chunkerkernevalmat(chnkr + ii*[d;0], fkern_0s, chnkr, struct('forcesmooth',true)), 1, chnkr.npt,   chnkr.npt);
end

M(:,2:2:end,1:2:end) = M(:,2:2:end,1:2:end) + M2;
sys = M3 + M;

% solve
ising = 1;
src    = struct('r', [0; -1.5]);
bskern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_bcs', kappa, d, sn, s0_l, sn_l, l, ising, nu);
skern  = @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',                   kappa, d, sn, s0_l, sn_l, l, ising);

rhs = -bskern(src, chnkr);
sol = squeeze(sys) \ rhs;

% check: total field vanishes exterior => uscat + uin = 0
targ    = struct('r', [3*d/4; 1.5]);
ikern   = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval', kappa, d, sn, s0_l, sn_l, l, 0, nu);
ikern_0 = @(s,t) chnk.flex2d.kern(    zk, s, t, 'supported_plate_eval', nu);

wts   = repmat(chnkr.wts(:).', 2, 1);
uscat = (ikern(chnkr, targ) .* wts(:).' + chunkerkernevalmat(chnkr, ikern_0, targ, opts_l)) * sol;
uin   = skern(src, targ);
assert(abs(uscat + uin) < 1e-8)

end


function [r,d,d2] = sin_func(t,d,A)
omega = 2*pi/d;
r = [t, A*sin(omega*t)].';
d = [ones(length(t),1), omega*A*cos(omega*t)].';
d2 = [zeros(length(t),1), -omega^2*A*sin(omega*t)].';
end