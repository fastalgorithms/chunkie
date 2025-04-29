quasiperiodicTest0();
quasiperiodicTest1();

function quasiperiodicTest0()
% test the representation

% problem parameters
d= 1;
zk = 1;
kappa = pi-0.1-1i;
coefs = [1;0.5];
coefa = [1,0.3;-1i,0.2];

src = []; src.r = [0;-1]; src.n = [1;-2];
targ = []; targ.r = [1.1;0.3]; targ.n = [-1;0.3];

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
assert(abs(svalshift - exp(1i*kappa*d)*sval) <  1e-12)

% test combined
assert(abs((coefs(1)*dval + coefs(2)*sval)-cval) <  1e-12)
assert(abs((coefs(1)*dpval + coefs(2)*spval)-cpval) <  1e-12)

% test all
assert(norm([coefa(1,1)*dval, coefa(1,2)*sval; coefa(2,1)*dpval, ...
    coefa(2,2)*spval] - allval) <  1e-12)

% test transmission representiatoon
assert(norm([coefs(1)*dval , coefs(2)*sval]-trval) <  1e-12)
assert(norm([coefs(1)*dpval , coefs(2)*spval]-trpval) <  1e-12)

% test gradients
assert(norm((coefs(1)*dgval + coefs(2)*sgval)-cgval) <  1e-12)
assert(norm([coefs(1)*dgval , coefs(2)*sgval]-trgval) <  1e-12)

% test c2trans
assert(norm([coefs(1)*dval + coefs(2)*sval; coefs(1)*dpval + coefs(2)*spval]-c2trval) <  1e-12)
end


function quasiperiodicTest1()
% test an integral equation


% problem parameters
d= 8;
zk = 1;
kappa = .05+0.1i;

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


% ploting
Lplot = d/1.5;
nplot = 80;
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
% maxu = .6;
figure(1);clf
subplot(1,3,1)
h = pcolor(X,Y,reshape(imag(uin),nplot,[])); set(h,'edgecolor','none')
title('$u_{\rm{in}}$','Interpreter','latex'); set(gca,'fontsize',14)
colorbar; clim([-maxu,maxu])
hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off

subplot(1,3,2)
h= pcolor(X,Y,reshape(imag(uscat),nplot,[])); set(h,'edgecolor','none')
title('$u_{\rm{scat}}$','Interpreter','latex'); set(gca,'fontsize',14)
colorbar; clim([-maxu,maxu])
hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off

subplot(1,3,3)
h = pcolor(X,Y,reshape(log10(abs(utot)),nplot,[])); set(h,'edgecolor','none')
title('$\log_{10} $ error','Interpreter','latex'); set(gca,'fontsize',14)
colorbar; % clim([-maxu,maxu])
hold on, plot(chnkr,'.'), plot(src.r(1,:),src.r(2,:),'o','LineWidth',2), hold off


end

function [r,d,d2] = sin_func(t,d,A)
omega = 2*pi/d;
r = [t, A*sin(omega*t)].';
d = [ones(length(t),1), omega*A*cos(omega*t)].';
d2 = [zeros(length(t),1), -omega^2*A*sin(omega*t)].';
end