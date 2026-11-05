quasiperiodicTest0();
quasiperiodicTest1();
quasiperiodicTest2();
quasiperiodicTest_flex_clamped();
quasiperiodicTest_flex_free();
quasiperiodicTest_flex_supported();
quasiperiodicTest_lap_algebraic();
quasiperiodicTest_flex_algebraic();
quasiperiodicTest_constructors();
quasiperiodicTest_reference();

function quasiperiodicTest0()
% algebraic consistency checks for helm2dquas kernels (Bloch shift, Green's function identity, ising flag)

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

allkern = kernel('hq','all',zk,kappa,d,coefa);
trkern = kernel('hq','trans_rep',zk,kappa,d,coefs);
trpkern = kernel('hq','trans_rep_prime',zk,kappa,d,coefs);

% gradient kernels have no kernel(...) constructor
sgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'sgrad',quas_param));
dgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'dgrad',quas_param));
cgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'cgrad',quas_param,coefs));
trgkern = kernel(@(s,t) chnk.helm2dquas.kern(zk,s,t,'trans_rep_grad',quas_param,coefs));

c2trkern = kernel('hq','c2trans',zk,kappa,d,coefs);

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

% test dp directly from the hessian of the quasiperiodic Green's function:
%   D'(x,y) = -hess . (n_src (x) n_targ)
sn_q = quas_param.sn; l_q = quas_param.l;
[~,~,hess_dp] = chnk.helm2dquas.green(src.r, targ.r, zk, kappa, d, sn_q, l_q, 1);
nxs = repmat(src.n(1,:), nkappa*nt, 1);
nys = repmat(src.n(2,:), nkappa*nt, 1);
nxt = repmat(reshape(targ.n(1,:),1,nt),nkappa,1); nxt = repmat(nxt(:),1,ns);
nyt = repmat(reshape(targ.n(2,:),1,nt),nkappa,1); nyt = repmat(nyt(:),1,ns);
dp_ref = -(hess_dp(:,:,1).*nxs.*nxt + hess_dp(:,:,2).*(nys.*nxt + nxs.*nyt) ...
    + hess_dp(:,:,3).*nys.*nyt);
assert(norm(dp_ref - dpval) < 1e-11)

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

% test nsub reconstruction for s/sp/dp (value/grad/hess branches of green):
%   kern(ising=0,nsub) + sum_{ii=-nsub}^{nsub} alpha^ii K_fs(targ-src-ii*[d;0])
%   = kern(ising=1,nsub=0).
% nsub must be <= l (default l=2 in the 'hq' factory).
nsub_t = 2;
alpha = exp(1i*kappa(:)*d);
for kt = {'s','sp','dp'}
    type = kt{1};
    kquas = kernel('hq',type,zk,kappa,d,[],struct('nsub',nsub_t),0);
    kfs   = kernel('h',type,zk);
    recon = kquas.eval(src,targ);
    for ii = -nsub_t:nsub_t
        src_ii = src; src_ii.r = src.r + ii*[d;0];
        Kfs = kfs.eval(src_ii,targ);   % free-space, (nt,ns)
        Kfs = repmat(reshape(Kfs,1,nt,ns), nkappa,1,1);
        Kfs = reshape(Kfs, nkappa*nt, ns);
        recon = recon + repmat(alpha.^ii, nt, 1).*Kfs;
    end
    kfull = kernel('hq',type,zk,kappa,d).eval(src,targ);
    assert(norm(recon - kfull) < 1e-11)
end

allkern1a = kernel('hq','all',zk,kappa,d,coefa,[],0);
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

% plotting (disabled in regression test; uncomment to visualize)
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
% toc
% targup = targ(:,iup.');
%
% tic;
% uin = NaN*zeros(ntarg,1)+NaN*1i;
% src.n = [0;1];
% uin(iup) = skern.eval(src,struct("r",targup));
% uscat = NaN*zeros(ntarg,1)+NaN*1i;
% opts = []; opts.forcesmooth = false; opts.eps = 1e-6;
% opts.forcepquad = true; opts.side = 'e';
% uscat(iup) = chunkerkerneval(chnkr,dkern,soln,targup,opts);
% toc
% %
% utot = uin-uscat;
% maxu = max(abs(uin));
% figure(1);clf
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

end

function quasiperiodicTest_flex_clamped()
% test an integral equation for the quasiperiodic clamped plate problem

% problem parameters
d = 2; zk = 7; kappa = 0.2-0.1*1i;

% setup geometry
nch = 2^3; A = .2;
cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
chnkr = reverse(chunkerfuncuni(@(t) sin_func(t,d,A), nch, cparams));
curv = signed_curvature(chnkr); curv = curv(:);

opts_s = []; opts_s.sing = 'smooth'; opts_s.quad = 'native';
opts_l = []; opts_l.sing = 'log';

% build system: quasi periodic part uses smooth quadrature (ising=0): add
% the self-interaction and nsub nearest periodic copies to ising=0 kernel
% and compare with original value
ising = 0; nsub = 1;
alpha = exp(1i*kappa*d);

fkern   = kernel('fq', 'clamped_plate', zk, kappa, d, [], struct('nsub',nsub), ising);
fkern_0 = @(s,t) chnk.flex2d.kern(    zk, s, t, 'clamped_plate');

sys   = reshape(chunkermat(chnkr, fkern,   opts_s), 1, 2*chnkr.npt, 2*chnkr.npt);
sys_0 = reshape(chunkermat(chnkr, fkern_0, opts_l), 1, 2*chnkr.npt, 2*chnkr.npt);
for ii = [-nsub:-1, 1:nsub]
    sys_0 = sys_0 + alpha^ii .* reshape(chunkerkernevalmat(chnkr + ii*[d;0], fkern_0, chnkr), 1, 2*chnkr.npt, 2*chnkr.npt);
end

sys = sys + sys_0 - reshape(0.5*eye(2*chnkr.npt), 1, 2*chnkr.npt, 2*chnkr.npt);
sys(:,2:2:end,1:2:end) = sys(:,2:2:end,1:2:end) + reshape(curv .* eye(chnkr.npt), 1, chnkr.npt, chnkr.npt);

% solve
src    = struct('r', [0; -1.5]);
bskern = kernel('fq', 'clamped_plate_bcs', zk, kappa, d);
skern  = kernel('fq', 's',                 zk, kappa, d);

rhs = -bskern.eval(src, chnkr);
sol = squeeze(sys) \ rhs;

% check: total field vanishes exterior => uscat + uin = 0
targ    = struct('r', [3*d/4; 1.5]);
ikern   = kernel('fq', 'clamped_plate_eval', zk, kappa, d, [], [], 0);
ikern_0 = @(s,t) chnk.flex2d.kern(    zk, s, t, 'clamped_plate_eval');

wts    = repmat(chnkr.wts(:).', 2, 1);
uscat  = (ikern.eval(chnkr, targ) .* wts(:).' + chunkerkernevalmat(chnkr, ikern_0, targ, opts_l)) * sol;
uin    = skern.eval(src, targ);
assert(abs(uscat + uin) < 1e-8)

end


function quasiperiodicTest_flex_free()
% test an integral equation for the quasiperiodic free plate problem

% problem parameters
d = 2; zk = 7; kappa = 0.2-0.1*1i; nu = 0.3;

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

qopts   = struct('nsub', nsub);
fkern   = kernel('fq', 'free_plate', zk, kappa, d, nu, qopts, ising);
double  = kernel('lq', 'd',    kappa, d, qopts, ising);
hilbert = kernel('lq', 'hilb', kappa, d, qopts, ising);
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
src    = struct('r', [0; -1.5]);
bskern = kernel('fq', 'free_plate_bcs', zk, kappa, d, nu);
skern  = kernel('fq', 's',              zk, kappa, d);

rhs = -bskern.eval(src, chnkr);
sol = squeeze(sys) \ rhs;

% build combined density [phi; H*phi; psi] for free_plate_eval
dens_comb = zeros(3*chnkr.npt, size(rhs,2));
dens_comb(1:3:end,:) = sol(1:2:end,:);
dens_comb(2:3:end,:) = squeeze(H) * sol(1:2:end,:);
dens_comb(3:3:end,:) = sol(2:2:end,:);

% check: total field vanishes exterior => uscat + uin = 0
targ    = struct('r', [3*d/4; 1.5]);
ikern   = kernel('fq', 'free_plate_eval', zk, kappa, d, nu, [], 0);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);

wts   = repmat(chnkr.wts(:).', 3, 1);
uscat = (ikern.eval(chnkr, targ) .* wts(:).' + chunkerkernevalmat(chnkr, ikern_0, targ, opts_l)) * dens_comb;
uin   = skern.eval(src, targ);
assert(abs(uscat + uin) < 1e-8)

end


function quasiperiodicTest_flex_supported()
% test an integral equation for the quasiperiodic supported plate problem

% problem parameters
d = 2; zk = 7; kappa = 0.2-0.1*1i; nu = 0.3;

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

% build system: quasi periodic part uses smooth quadrature (ising=0): add
% the self-interaction and nsub nearest periodic copies to ising=0 kernel
% and compare with original value
ising = 0; nsub = 1;
alpha = exp(1i*kappa*d);
c0    = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

fkern    = kernel('fq', 'supported_plate', zk, kappa, d, nu, struct('nsub',nsub), ising);
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
src    = struct('r', [0; -1.5]);
bskern = kernel('fq', 'supported_plate_bcs', zk, kappa, d, nu);
skern  = kernel('fq', 's',                   zk, kappa, d);

rhs = -bskern.eval(src, chnkr);
sol = squeeze(sys) \ rhs;

% check: total field vanishes exterior => uscat + uin = 0
targ    = struct('r', [3*d/4; 1.5]);
ikern   = kernel('fq', 'supported_plate_eval', zk, kappa, d, nu, [], 0);
ikern_0 = @(s,t) chnk.flex2d.kern(    zk, s, t, 'supported_plate_eval', nu);

wts   = repmat(chnkr.wts(:).', 2, 1);
uscat = (ikern.eval(chnkr, targ) .* wts(:).' + chunkerkernevalmat(chnkr, ikern_0, targ, opts_l)) * sol;
uin   = skern.eval(src, targ);
assert(abs(uscat + uin) < 1e-8)

end


function quasiperiodicTest_lap_algebraic()
% algebraic consistency checks for lap2dquas kernels (Bloch shift, Green's function identity, ising flag, nsub)

% problem parameters
d = 1;
kappa = [pi-0.1-1i, 0.2, -0.1, 0.25i];
nkappa = length(kappa);

src = []; src.r = [[0;-1.1],[1;-1],[0.1;-0.3]]; src.n = [[1;-2],[2;-1],[1;0]];
targ = []; targ.r = [[1.1;0.3],[2;0]]; targ.n = [[-1;0.3],[0.1;1]];
ns = size(src.r,2);
nt = size(targ.r,2);

% kernels via @kernel
skern  = kernel('lq', 's',  kappa, d);
spkern = kernel('lq', 'sp', kappa, d);
dkern  = kernel('lq', 'd',  kappa, d);
dpkern = kernel('lq', 'dp', kappa, d);
stkern = kernel('lq', 'st', kappa, d);

quas_param = skern.params.quas_param;
s0q = quas_param.s0; snq = quas_param.sn; lq = quas_param.l;

hilbkern  = kernel('lq', 'hilb',      kappa, d);
hilbpkern = kernel('lq', 'hilbprime', kappa, d);

% evaluate all kernels
sval  = skern.eval(src,targ);
spval = spkern.eval(src,targ);
dval  = dkern.eval(src,targ);
dpval = dpkern.eval(src,targ);
stval = stkern.eval(src,targ);

hilbval  = hilbkern.eval(src,targ);
hilbpval = hilbpkern.eval(src,targ);

% test quasi-periodicity: shift target by d => multiply by exp(i*kappa*d)
targ_s = targ; targ_s.r = targ.r + [d;0];
svalshift = skern.eval(src,targ_s);
assert(norm(svalshift - repmat(exp(1i*kappa(:)*d),nt,1).*sval) < 1e-12)

% test hilb directly from definition:
%   hilb(x,y) = 2 * grad_x G(x,y) . tau_y  where tau_y = (n_y^src, -n_x^src)
[~,grad] = chnk.lap2dquas.green(src.r, targ.r, kappa, d, s0q, snq, lq, 1);
% grad is (nkappa*nt, ns, 2); src.n gives tau_src = (n_y, -n_x)
tau_x = repmat(src.n(2,:), nkappa*nt, 1);   % n_y component
tau_y = repmat(-src.n(1,:), nkappa*nt, 1);  % -n_x component
hilb_ref = 2*(grad(:,:,1).*tau_x + grad(:,:,2).*tau_y);
assert(norm(hilb_ref - hilbval) < 1e-11)

% test st (tangential derivative of single layer) directly from definition:
%   st(x,y) = grad_x G . tau_x  where tau_x = (-n_y^targ, n_x^targ)
nxt = repmat(reshape(targ.n(1,:),1,nt),nkappa,1); nxt = nxt(:);
nyt = repmat(reshape(targ.n(2,:),1,nt),nkappa,1); nyt = nyt(:);
nxt = repmat(nxt,1,ns); nyt = repmat(nyt,1,ns);
st_ref = -grad(:,:,1).*nyt + grad(:,:,2).*nxt;
assert(norm(st_ref - stval) < 1e-11)

% test hilbprime (tangential derivative of Hilbert transform) directly from
% definition, using the hessian of the quasiperiodic Green's function.
[~,~,hess_h] = chnk.lap2dquas.green(src.r, targ.r, kappa, d, s0q, snq, lq, 1);
nx_s = repmat(src.n(1,:), nkappa*nt, 1);
ny_s = repmat(src.n(2,:), nkappa*nt, 1);
hilbp_ref = 2*(-(hess_h(:,:,1).*ny_s - hess_h(:,:,2).*nx_s).*nyt ...
    + (hess_h(:,:,2).*ny_s - hess_h(:,:,3).*nx_s).*nxt);
assert(norm(hilbp_ref - hilbpval) < 1e-11)

% test ising flag against the Green's function definition directly.
% S kernel is G; D kernel is -dG/dn_src = -(G_grad . n_src). The ising=1
% values are already covered by the Bloch-shift and hilb checks above
% (sval == green(...,1)), so here we pin down the ising=0 kernels.
sval0 = kernel('lq', 's', kappa, d, [], 0).eval(src,targ);
dval0 = kernel('lq', 'd', kappa, d, [], 0).eval(src,targ);

[G0, G0_grad] = chnk.lap2dquas.green(src.r, targ.r, kappa, d, s0q, snq, lq, 0);
nx_fs = repmat(src.n(1,:), nkappa*nt, 1);
ny_fs = repmat(src.n(2,:), nkappa*nt, 1);
dG0 = -(G0_grad(:,:,1).*nx_fs + G0_grad(:,:,2).*ny_fs);

assert(norm(sval0 - G0) < 1e-12)
assert(norm(dval0 - dG0) < 1e-12)

% ising=1 D kernel against -dG/dn_src, reusing grad = green(...,1) from above
dG1 = -(grad(:,:,1).*nx_fs + grad(:,:,2).*ny_fs);
assert(norm(dval - dG1) < 1e-12)

% test nsub reconstruction for s/sp/dp (value/grad/hess branches of green):
%   kern(ising=0,nsub) + sum_{ii=-nsub}^{nsub} alpha^ii K_fs(targ-src-ii*[d;0])
%   = kern(ising=1,nsub=0).
% The kernel('lq',...) wrapper also exercises nsub threading via quad_opts.
% nsub must be <= l (here l = 2).
nsub_t = 2;
alpha = exp(1i*kappa(:)*d);
for kt = {'s','sp','dp'}
    type = kt{1};
    kquas = kernel('lq',type,kappa,d,struct('nsub',nsub_t),0);
    kfs   = kernel('l',type);
    recon = kquas.eval(src,targ);
    for ii = -nsub_t:nsub_t
        src_ii = src; src_ii.r = src.r + ii*[d;0];
        Kfs = kfs.eval(src_ii,targ);   % free-space, (nt,ns)
        Kfs = repmat(reshape(Kfs,1,nt,ns), nkappa,1,1);
        Kfs = reshape(Kfs, nkappa*nt, ns);
        recon = recon + repmat(alpha.^ii, nt, 1).*Kfs;
    end
    kfull = kernel('lq',type,kappa,d).eval(src,targ);
    assert(norm(recon - kfull) < 1e-11)
end

end


function quasiperiodicTest_flex_algebraic()
% algebraic consistency checks for flex2dquas kernels (Bloch shift, Green's function identity, ising flag, nsub)

% problem parameters
d = 1; zk = 3;
kappa = [pi-0.1-1i, 0.2, -0.1, 0.25i];
nkappa = length(kappa);

src = []; src.r = [[0;-1.1],[1;-1],[0.1;-0.3]];
src.n = [[1;-2],[2;-1],[1;0]]; src.n = src.n ./ vecnorm(src.n);
src.d = [[2;1],[1;2],[0;1]];
targ = []; targ.r = [[1.1;0.3],[2;0]];
targ.n = [[-1;0.3],[0.1;1]]; targ.n = targ.n ./ vecnorm(targ.n);
targ.d = [[1;-1],[1;1]]; targ.d2 = [[0;1],[0;-1]];
ns = size(src.r,2);
nt = size(targ.r,2);

% kernels via @kernel
K = @(type) kernel('fq', type, zk, kappa, d);
Sn = K('s').params.quas_param.sn; l = K('s').params.quas_param.l;

sval  = K('s').eval(src,targ);
dval  = K('d').eval(src,targ);
spval = K('sp').eval(src,targ);

% test quasi-periodicity: shift target by d => multiply by exp(i*kappa*d).
% Check s (value), d and sp (grad branch) so the gradient phase is exercised.
targ_s = targ; targ_s.r = targ.r + [d;0];
bloch = repmat(exp(1i*kappa(:)*d),nt,1);
assert(norm(K('s').eval(src,targ_s)  - bloch.*sval)  < 1e-12)
assert(norm(K('d').eval(src,targ_s)  - bloch.*dval)  < 1e-12)
assert(norm(K('sp').eval(src,targ_s) - bloch.*spval) < 1e-12)

% test that flex S kernel = G_flex / (2*zk^2)
% where G_flex is the raw quasiperiodic flexural Green's function
Gflex = chnk.flex2dquas.green(src.r, targ.r, zk, kappa, d, Sn, l, 1);
sval_ref = reshape(Gflex, nkappa*nt, ns) / (2*zk^2);
assert(norm(sval - sval_ref) < 1e-11)

% test that G_flex = G_helm(zk) - G_helm(i*zk)
% (the defining relation of the flexural Green's function)
shkern_zk  = kernel('hq', 's', zk,    kappa, d);
shkern_izk = kernel('hq', 's', 1i*zk, kappa, d);
Gh_zk  = chnk.helm2dquas.green(src.r, targ.r, zk,    kappa, d, shkern_zk.params.quas_param.sn,  shkern_zk.params.quas_param.l,  1);
Gh_izk = chnk.helm2dquas.green(src.r, targ.r, 1i*zk, kappa, d, shkern_izk.params.quas_param.sn, shkern_izk.params.quas_param.l, 1);
assert(norm(Gflex - (Gh_zk - Gh_izk)) < 1e-11)

% test ising flag: ising=0 (periodic images) + free-space singular part = ising=1
K0  = @(type) kernel('fq', type, zk, kappa, d, [], [], 0);
K1  = K;

sval0  = K0('s').eval(src,targ);
dval0  = K0('d').eval(src,targ);
spval0 = K0('sp').eval(src,targ);

% free-space flex kern (chnk.flex2d.kern)
sval_fs  = kernel(@(s,t) chnk.flex2d.kern(zk,s,t,'s')).eval(src,targ);
dval_fs  = kernel(@(s,t) chnk.flex2d.kern(zk,s,t,'d')).eval(src,targ);
spval_fs = kernel(@(s,t) chnk.flex2d.kern(zk,s,t,'sp')).eval(src,targ);

sval_fs  = repmat(reshape(sval_fs, 1,nt,ns), nkappa,1,1);
dval_fs  = repmat(reshape(dval_fs, 1,nt,ns), nkappa,1,1);
spval_fs = repmat(reshape(spval_fs,1,nt,ns), nkappa,1,1);

assert(norm(sval0  + reshape(sval_fs, nkappa*nt,ns)  - K1('s').eval(src,targ))  < 1e-11)
assert(norm(dval0  + reshape(dval_fs, nkappa*nt,ns)  - K1('d').eval(src,targ))  < 1e-11)
assert(norm(spval0 + reshape(spval_fs,nkappa*nt,ns)  - K1('sp').eval(src,targ)) < 1e-11)

% test nsub reconstruction for s/sp/d (value/grad/hess branches of green;
% dp is not a standalone flex kernel). Goes through the kernel('fq',...)
% wrapper, so this also value-tests nsub threading via quad_opts:
%   kern(ising=0,nsub) + sum_{ii=-nsub}^{nsub} alpha^ii K_fs(targ-src-ii*[d;0])
%   = kern(ising=1,nsub=0).
% nsub must be <= l (here l = 2).
nsub_t = 2;
alpha = exp(1i*kappa(:)*d);
for kt = {'s','sp','d'}
    type = kt{1};
    Kfs_h = kernel(@(s,t) chnk.flex2d.kern(zk,s,t,type));
    recon = kernel('fq',type,zk,kappa,d,[],struct('nsub',nsub_t),0).eval(src,targ);
    for ii = -nsub_t:nsub_t
        src_ii = src; src_ii.r = src.r + ii*[d;0];
        Kfs = Kfs_h.eval(src_ii, targ);   % free-space, (nt,ns)
        Kfs = repmat(reshape(Kfs,1,nt,ns), nkappa,1,1);
        Kfs = reshape(Kfs, nkappa*nt, ns);
        recon = recon + repmat(alpha.^ii, nt, 1).*Kfs;
    end
    kfull = K(type).eval(src,targ);
    assert(norm(recon - kfull) < 1e-11)
end

end


function [r,d,d2] = sin_func(t,d,A)
omega = 2*pi/d;
r = [t, A*sin(omega*t)].';
d = [ones(length(t),1), omega*A*cos(omega*t)].';
d2 = [zeros(length(t),1), -omega^2*A*sin(omega*t)].';
end


function quasiperiodicTest_constructors()
% regression tests that go through the kernel(...) constructors and compare
% against independent references (spectral dual sum, brute-force images)

d = 1; zk = 1.3; kappa = 0.3; alpha = exp(1i*kappa*d);
src = []; src.r = [0.1;-0.2]; src.n = [0.6;0.8];
src.d = [0.8;-0.6]; src.d2 = [0.1;0.2];
targ = []; targ.r = [0.5;0.3]; targ.n = [-0.8;0.6];
targ.d = [0.6;0.8]; targ.d2 = [-0.1;0.3];
rx = targ.r(1)-src.r(1); ry = targ.r(2)-src.r(2);

% Laplace lattice coefficients stored by flex2dquas match l (used by the
% free plate hilbert transforms); check hilb against the dual sum gradient
l = 2;
fk = kernel('fq','free_plate',zk,kappa,d,[0.3 0]);
lk = kernel('lq','hilb',kappa,d);
assert(norm(fk.params.quas_param.sn_l - lk.params.quas_param.sn) < 1e-14);
assert(norm(fk.params.quas_param.s0_l - lk.params.quas_param.s0) < 1e-14);
[~,g] = chnk.lap2dquas.quasi_dual_sum(rx,ry,0,kappa,d);
hilbref = 2*(g(1,1,1)*src.n(2) - g(1,1,2)*src.n(1));
hilb = chnk.lap2dquas.kern(src,targ,'hilb',kappa,d,fk.params.quas_param.s0_l, ...
    fk.params.quas_param.sn_l,l,1);
assert(abs(hilb-hilbref) < 1e-12);
assert(abs(lk.eval(src,targ)-hilbref) < 1e-12);

% ising = 0 kernels are finite and continuous at coincident points
for type = {'s','d','sp','dp'}
    kh = kernel('hq',type{1},zk,kappa,d,[],[],0);
    v = kh.eval(src,src);
    assert(all(isfinite(v(:))));
    s2 = src; s2.r = src.r + 1e-7*[1;1];
    assert(abs(v - kh.eval(s2,src)) < 1e-5);
end

% ising = 0 kernels use smooth quadrature, and all carry quas_param
for type = {'s','d','sp','dp','st','hilb','hilbprime'}
    k = kernel('lq',type{1},kappa,d,[],0);
    assert(strcmp(k.sing,'smooth') && isfield(k.params,'quas_param'));
end
for type = {'s','d','sp','clamped_plate','free_plate','supported_plate'}
    k = kernel('fq',type{1},zk,kappa,d,0.3,[],0);
    assert(strcmp(k.sing,'smooth') && isfield(k.params,'quas_param'));
    assert(k.params.quas_param.kappa == kappa);
end

% combined kernels without explicit coefs fall back to defaults
ck = kernel('hq','c',zk,kappa,d);
assert(abs(ck.eval(src,targ) - ...
    kernel('hq','c',zk,kappa,d,[1,1i]).eval(src,targ)) < 1e-14);
kernel('hq','cp',zk,kappa,d);

% block kernels respect ising and have opdims consistent with eval
kap2 = [0.3 0.4]; coefa = [1,0.3;-1i,0.2]; coefs = [1;0.5];
s0k = kernel('hq','s',zk,kap2,d,[],[],0).eval(src,targ);
d0k = kernel('hq','d',zk,kap2,d,[],[],0).eval(src,targ);
types = {'all','trep','trep_p','c2t'};
cfs = {coefa,coefs,coefs,coefa};
for it = 1:numel(types)
    k = kernel('hq',types{it},zk,kap2,d,cfs{it},[],0);
    v = k.eval(src,targ);
    assert(isequal(size(v), k.opdims));
    assert(strcmp(k.sing,'smooth'));
end
va = kernel('hq','all',zk,kap2,d,coefa,[],0).eval(src,targ);
va = reshape(va,numel(kap2),2,[],2); % rows ordered (kappa, component, target)
assert(norm(reshape(va(:,1,:,2),[],1) - coefa(1,2)*s0k(:)) < 1e-14);
vt = kernel('hq','trep',zk,kap2,d,coefs,[],0).eval(src,targ);
assert(norm(vt(:,1) - coefs(1)*d0k) < 1e-14);

% clamped plate: default constructor = periodic part (nsub=0) + free space
ck = kernel('fq','clamped_plate',zk,kappa,d);
Kp = chnk.flex2dquas.kern(zk,src,targ,'clamped_plate',kappa,d, ...
    ck.params.quas_param.sn,[],[],l,0,0);
Kf = chnk.flex2d.kern(zk,src,targ,'clamped_plate');
assert(norm(ck.eval(src,targ) - (Kp+Kf)) < 1e-12*norm(Kp+Kf));

% quad_opts.nsub reaches the plate kernels
qo = []; qo.nsub = 1;
ck1 = kernel('fq','clamped_plate',zk,kappa,d,[],qo,0);
Kp1 = chnk.flex2dquas.kern(zk,src,targ,'clamped_plate',kappa,d, ...
    ck.params.quas_param.sn,[],[],l,0,1);
assert(norm(ck1.eval(src,targ) - Kp1) < 1e-14);
fk1 = kernel('fq','free_plate',zk,kappa,d,[0.3 0],qo,0);
Fp1 = chnk.flex2dquas.kern(zk,src,targ,'free_plate',kappa,d, ...
    fk1.params.quas_param.sn,fk1.params.quas_param.s0_l,fk1.params.quas_param.sn_l,l,0,0.3,1);
assert(norm(fk1.eval(src,targ) - Fp1) < 1e-14);

% 'st' is accepted by the low-level Laplace kern
[s0,sn] = chnk.lap2dquas.latticecoefs(1:40,d,kappa,l+1);
assert(abs(chnk.lap2dquas.kern(targ,src,'st',kappa,d,s0,sn,l,1) - ...
    kernel('lq','st',kappa,d).eval(targ,src)) < 1e-14);

end


function quasiperiodicTest_reference()
% compare against references independent of the lattice sums:
%  - green (all derivative orders) vs the spectral sum
%  - helm kernels vs explicit image sums (complex zk, so images decay)
%  - flex kernels: ising=1 minus ising=0 (nsub=1) vs the free-space
%    kernels at the 3 nearest images

d = 1; kaps = [0.3, -1.1+0.2i];
src = []; src.r = [0.1;0.2];
rx = [0.37, -0.49, 1.8, -2.45]; ry = [0.42, -0.3, -0.9, 0.15];
targ = []; targ.r = src.r + [rx; ry];

zk = 1.3;
kl = kernel('lq','s',kaps,d).params.quas_param;
kh = kernel('hq','s',zk,kaps,d).params.quas_param;
kf = kernel('fq','s',zk,kaps,d).params.quas_param;
gl = cell(1,3); gh = cell(1,3); gf = cell(1,7);
[gl{:}] = chnk.lap2dquas.green(src.r,targ.r,kaps,d,kl.s0,kl.sn,kl.l,1);
[gh{:}] = chnk.helm2dquas.green(src.r,targ.r,zk,kaps,d,kh.sn,kh.l,1);
[gf{:}] = chnk.flex2dquas.green(src.r,targ.r,zk,kaps,d,kf.sn,kf.l,1);
for o = 0:6
    rf = spectral(rx,ry,zk,kaps,d,o) - spectral(rx,ry,1i*zk,kaps,d,o);
    assert(relerr(gf{o+1},rf) < 1e-10);
    if o <= 2
        assert(relerr(gl{o+1},spectral(rx,ry,0,kaps,d,o)) < 1e-12);
        assert(relerr(gh{o+1},spectral(rx,ry,zk,kaps,d,o)) < 1e-12);
    end
end

% kernel types, through the constructors
src = []; src.r = [[0.1;0.2],[-0.3;-0.1]]; src.n = [[0.6;0.8],[1;0]];
src.d = [[0.8;-0.6],[0;-2]]; src.d2 = [[0.1;0.2],[0.3;-0.2]];
targ = []; targ.r = [[0.5;0.6],[-1.6;-0.4],[2.3;0.1]];
targ.n = [[-0.8;0.6],[0;1],[0.6;-0.8]];
targ.d = [[0.6;0.8],[2;0],[-1.6;-1.2]];
targ.d2 = [[-0.1;0.3],[0.2;0.1],[0;0.4]];
src.data = [[0.1;0.3],[-0.2;0.4]]; targ.data = [[0.2;0.1],[0.3;-0.2],[0.1;0.5]];
kap = 0.3; alpha = exp(1i*kap*d); nu = 0.3;
coefs = [1;0.5]; coefa = [1,0.3;-1i,0.2];
qo = []; qo.nsub = 1;

zkc = 2+2i;
helm = {'s',[],'d',[],'sp',[],'dp',[],'c',coefs,'cp',coefs,'all',coefa, ...
    'trep',coefs,'trep_p',coefs,'c2t',coefa};
for it = 1:2:numel(helm)
    ty = helm{it}; cf = helm{it+1};
    kfs = @(s,t) chnk.helm2d.kern(zkc,s,t,ty,cf);
    Q1 = kernel('hq',ty,zkc,kap,d,cf).eval(src,targ);
    Q0 = kernel('hq',ty,zkc,kap,d,cf,qo,0).eval(src,targ);
    assert(relerr(Q1,imagesum(kfs,src,targ,-20:20,alpha,d)) < 1e-12);
    assert(relerr(Q0,imagesum(kfs,src,targ,[-20:-2,2:20],alpha,d)) < 1e-12);
end
qp = kernel('hq','s',zkc,kap,d).params.quas_param;
for ty = {'sgrad','dgrad','cgrad','trep_g'}
    Q = chnk.helm2dquas.kern(zkc,src,targ,ty{1},qp,coefs);
    R = imagesum(@(s,t) chnk.helm2d.kern(zkc,s,t,ty{1},coefs), ...
        src,targ,-20:20,alpha,d);
    assert(relerr(Q,R) < 1e-12);
end

flex = {'s',[],'d',[],'sp',[],'clamped_plate_bcs',[],'clamped_plate',[], ...
    'clamped_plate_eval',[],'free_plate_bcs',nu,'free_plate',nu, ...
    'free_plate_eval',nu,'supported_plate_bcs',nu,'supported_plate_eval',nu};
for it = 1:2:numel(flex)
    ty = flex{it}; cf = flex{it+1};
    if isempty(cf), args = {}; else, args = {cf}; end
    Q = kernel('fq',ty,zk,kap,d,cf,qo,1).eval(src,targ) ...
        - kernel('fq',ty,zk,kap,d,cf,qo,0).eval(src,targ);
    R = imagesum(@(s,t) chnk.flex2d.kern(zk,s,t,ty,args{:}), ...
        src,targ,-1:1,alpha,d);
    assert(relerr(Q,R) < 1e-12);
end

% supported plate: the free-space part is the log kernel plus the smooth
% kernel in the (2,1) block
Q = kernel('fq','supported_plate',zk,kap,d,nu,qo,1).eval(src,targ) ...
    - kernel('fq','supported_plate',zk,kap,d,nu,qo,0).eval(src,targ);
R = imagesum(@(s,t) chnk.flex2d.kern(zk,s,t,'supported_plate_log',nu), ...
    src,targ,-1:1,alpha,d);
R(2:2:end,1:2:end) = R(2:2:end,1:2:end) + imagesum(@(s,t) ...
    chnk.flex2d.kern(zk,s,t,'supported_plate_smooth',nu),src,targ,-1:1,alpha,d);
assert(relerr(Q,R) < 1e-12);

% transmission kernels: bcs_trx columns are the bcs applied to the source
% fields G, -G_x, G_xx + 2 G_yy, -(G_xxx + 2 G_xyy) (columns 4...1), and
% eval_trx rows are d_x^(k-1) of eval (checked by finite differences)
src.r = src.r + [0;-0.8];
h = 1e-3; e = [h;0];
sh = @(s,v) setfield(s,'r',s.r+v);
for pl = {{'clamped',[]},{'free',nu},{'supported',nu}}
    nm = pl{1}{1}; cf = pl{1}{2};
    B = kernel('fq',[nm '_plate_bcs'],zk,kap,d,cf);
    BT = kernel('fq',[nm '_plate_bcs_trx'],zk,kap,d,cf);
    bt = BT.eval(src,targ);
    btx = (BT.eval(sh(src,e),targ) - BT.eval(sh(src,-e),targ))/(2*h);
    bx = (B.eval(sh(src,e),targ) - B.eval(sh(src,-e),targ))/(2*h);
    assert(relerr(bt(:,4:4:end),B.eval(src,targ)) < 1e-14);
    assert(relerr(bt(:,3:4:end),bx) < 1e-3);
    assert(relerr(bt(:,1:4:end),btx(:,2:4:end)) < 1e-3);

    E = kernel('fq',[nm '_plate_eval'],zk,kap,d,cf);
    ET = kernel('fq',[nm '_plate_eval_trx'],zk,kap,d,cf);
    et = ET.eval(src,targ);
    etx = (ET.eval(src,sh(targ,e)) - ET.eval(src,sh(targ,-e)))/(2*h);
    assert(relerr(et(1:4:end,:),E.eval(src,targ)) < 1e-14);
    for k = 2:4
        assert(relerr(et(k:4:end,:),etx(k-1:4:end,:)) < 1e-3);
    end
end

% ising = 0 flex Green's function is finite at r = 0 to all orders
gf = cell(1,7);
[gf{:}] = chnk.flex2dquas.green([0;0],[0;0],zk,kap,d, ...
    kernel('fq','s',zk,kap,d).params.quas_param.sn,2,0);
assert(all(cellfun(@(g) all(isfinite(g(:))), gf)));

end

function e = relerr(a,b)
e = norm(a(:)-b(:))/norm(b(:));
end

function R = imagesum(kfs,src,targ,js,alpha,d)
% sum of free-space kernels over the source images js
R = 0;
for j = js
    s = src; s.r = src.r + [j*d;0];
    R = R + alpha^j*kfs(s,targ);
end
end

function D = spectral(rx,ry,zk,kappa,d,o)
% derivatives d^o/(dx^(o-j) dy^j), j = 0...o, of the quasiperiodic Green's
% function from its spectral (Fourier) series, size (nkappa, npt, o+1)
kappa = kappa(:); rx = rx(:).'; ry = ry(:).';
M = ceil((40 + 6*o)/min(abs(ry))*d/(2*pi)) + 5;
xi = kappa + 2*pi/d*reshape(-M:M,1,1,[]);
beta = sqrt(1i*(xi-zk)).*sqrt(-1i*(xi+zk));
f = exp(-beta.*abs(ry) + 1i*xi.*rx)./(2*beta)/d;
D = zeros(numel(kappa),numel(rx),o+1);
for j = 0:o
    D(:,:,j+1) = sum((1i*xi).^(o-j).*(-beta.*sign(ry)).^j.*f,3);
end
end
