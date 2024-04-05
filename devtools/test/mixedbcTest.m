% MIXEDBCTEST test the code with mixed boundary conditions, first test with
% mixed Dirichlet and Neumann BCs then test with Dirichlet and transmission
% conditions, which has variable opdims

clearvars; close all;

nverts = 4; 
vertsout = exp(1i*2*pi*(0:(nverts-1))/nverts);
vertsout = [real(vertsout);imag(vertsout)];

vertsin = .5*exp(1i*2*pi*(0:(nverts-1))/nverts+1i*pi/4);
vertsin = [real(vertsin);imag(vertsin)];

verts = [vertsout, vertsin];

iind = 1:nverts;
jind = 1:nverts;

iind = [iind iind];
jind = [jind jind + 1];
jind(jind>nverts) = 1;
svals = [-ones(1,nverts) ones(1,nverts)];
edge2verts = sparse(iind,jind,svals,nverts,nverts);
edge2verts = [edge2verts, 0*edge2verts; 0*edge2verts, edge2verts];

edir = 1:4; % indices of edges with Dirichlet conditions
eneu = 5:8; % indices of edges with Neumann conditions

fchnks = cell(1,size(edge2verts,1));

cparams = [];
cparams.nover = 2;
[cgrph] = chunkgraph(verts,edge2verts,fchnks,cparams);

vstruc = procverts(cgrph);
rgns = findregions(cgrph);
cgrph = balance(cgrph);

% dirichlet and neumann test
zk = 30;

kerns(length([edir,eneu]),length([edir,eneu])) = kernel();
kerns(edir,edir) = -2*kernel('helm', 'd', zk); % Dirichlet conditions
kerns(eneu,edir) = -2*kernel('helm', 'dp', zk); % Neumann conditions

kerns(edir,eneu) = -2*kernel('helm', 's', zk);
kerns(eneu,eneu) = -2*kernel('helm', 'sp', zk);


start = tic; sysmat = chunkermat(cgrph,kerns); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = eye(size(sysmat,1)) + sysmat;

fkernsrc = kernel('helm','s',zk);
sources = [1;1];
charges = [1];
srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = merge(cgrph.echnks(edir)).r(:,:); 
targinfo.d = merge(cgrph.echnks(edir)).d(:,:);
bdrydatad = fkernsrc.fmm(1e-12,srcinfo,targinfo,charges);

fkernsrc = kernel('helm','sp',zk);
targinfo = []; targinfo.r = merge(cgrph.echnks(eneu)).r(:,:); 
targinfo.d = merge(cgrph.echnks(eneu)).d(:,:);
targinfo.n = merge(cgrph.echnks(eneu)).n(:,:);
bdrydatan = fkernsrc.fmm(1e-12,srcinfo,targinfo,charges);
bdrydata = [bdrydatad;bdrydatan];

sol1 = sys\bdrydata;

cormat = chunkermat(cgrph,kerns,struct("corrections",true));
sysapply = @(sigma) bieapply(cgrph,kerns,1,sigma,cormat);

xapply1 = sys*bdrydata;
xapply2 = sysapply(bdrydata);
relerr = norm(xapply1-xapply2)/norm(xapply1);
fprintf('relative matrix free apply error %5.2e\n',relerr);
assert(relerr < 1e-10)


rmin = min(cgrph.r(:,:)'); rmax = max(cgrph.r(:,:)');
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 100;
xtarg = linspace(rmin(1)-0.1*xl,rmax(1)+0.1*xl,nplot); 
ytarg = linspace(rmin(2)-0.1*yl,rmax(2)+0.1*yl,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; in1 = chunkerinterior(merge(cgrph.echnks(edir)),{xtarg,ytarg});
in2 = chunkerinterior(merge(cgrph.echnks(eneu)),{xtarg,ytarg}); 
t1 = toc(start);
in = in1 & ~in2;
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential based on oversample boundary

start = tic;
fkernd = 2*kernel('helm', 'd', zk);
iddir = 1:merge(cgrph.echnks(edir)).npt;
uscat = chunkerkerneval(merge(cgrph.echnks(edir)),fkernd,sol1(iddir), ...
    targets(:,in));

fkern  = 2*kernel('helm', 's', zk);
idneu = (1:merge(cgrph.echnks(eneu)).npt) + merge(cgrph.echnks(edir)).npt;
uscat = uscat + chunkerkerneval(merge(cgrph.echnks(eneu)),fkern, ...
    sol1(idneu),targets(:,in)); 
t1 = toc(start);

fprintf('%5.2e s : time for kernel eval (for plotting)\n',t1)
fkernsrc = kernel('helm','s',zk);
targinfo = []; targinfo.r = targets(:,in); 
uin = fkernsrc.fmm(1e-12,srcinfo,targinfo,charges);
utot = uscat(:)+uin(:);

figure(1)
t = tiledlayout(1,2,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg));
zztarg(in) = uscat;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(cgrph,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg));
zztarg(in) = utot;
h=pcolor(xxtarg,yytarg,log10(abs((zztarg))));
set(h,'EdgeColor','none')
hold on
plot(cgrph,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
colorbar()

title('$\log10($error$)$','Interpreter','latex','FontSize',12)
relerr = max(abs(utot));
fprintf('relative field error %5.2e\n',relerr);
assert(relerr < 1e-4)

% dirichlet and transmission test
nregions = 2;
ks = [1.1;2.1]*30;


kerns(length([edir,eneu]),length([edir,eneu])) = kernel();
kerns(edir,edir) = -2*kernel('helm', 'd', ks(1));

kerntmp = @(s,t) -2*[chnk.helm2d.kern(ks(1),s,t,'d'); ...
    -chnk.helm2d.kern(ks(1),s,t,'dprime')];
kerns(eneu,edir) = (kerntmp);

kerntmp = @(s,t) -chnk.helm2d.kern(ks(1),s,t,'trans_rep',[-1,-1]);
kerns(edir,eneu) = (kerntmp);

cc1 = [-1,-1;1,1];
cc2 = cc1;
kerntmp = @(s,t) -(chnk.helm2d.kern(ks(1),s,t,'all',cc1)- ...
                 chnk.helm2d.kern(ks(2),s,t,'all',cc2));
kerns(eneu,eneu) = (kerntmp);


start = tic; sysmat = chunkermat(cgrph,kerns); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = eye(size(sysmat,1)) + sysmat;

fkernsrc = kernel('helm','s',ks(1));
sources = [.1;.1];
charges = [1];
srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = merge(cgrph.echnks(edir)).r(:,:); 
targinfo.d = merge(cgrph.echnks(edir)).d(:,:);
bdrydatad = fkernsrc.fmm(1e-12,srcinfo,targinfo,charges);

srcinfo.sources = sources;
srcinfo.charges = charges;
targinfo = []; 
targinfo.sources = merge(cgrph.echnks(eneu)).r(:,:); 
targinfo.n = merge(cgrph.echnks(eneu)).n(:,:);
U = hfmm2d(1e-12,ks(1),srcinfo,0,targinfo.sources,2);
bdrydatat = [U.pottarg; -sum(targinfo.n.*U.gradtarg,1)];
bdrydata = [bdrydatad;bdrydatat(:)];

sol1 = sys\bdrydata;

cormat = chunkermat(cgrph,kerns,struct("corrections",true));
sysapply = @(sigma) bieapply(cgrph,kerns,1,sigma,cormat);

xapply1 = sys*bdrydata;
xapply2 = sysapply(bdrydata);
relerr = norm(xapply1-xapply2)/norm(xapply1);
fprintf('relative matrix free apply error %5.2e\n',relerr);
assert(relerr < 1e-10)

start = tic;
fkernd = 2*kernel('helm', 'd', ks(1));
iddir = 1:merge(cgrph.echnks(edir)).npt;
uscat1 = chunkerkerneval(merge(cgrph.echnks(edir)),fkernd,sol1(iddir), ...
    targets(:,in));

kerntmp = @(s,t) -chnk.helm2d.kern(ks(1),s,t,'trans_rep',[1,1]);
fkern  = kernel(kerntmp);
idneu = (1:2*merge(cgrph.echnks(eneu)).npt) + merge(cgrph.echnks(edir)).npt;
uscat1 = uscat1 - chunkerkerneval(merge(cgrph.echnks(eneu)),fkern, ...
    sol1(idneu),targets(:,in));

kerntmp = @(s,t) -chnk.helm2d.kern(ks(2),s,t,'trans_rep',[1,1]);
fkern  = kernel(kerntmp);
uscat2 = chunkerkerneval(merge(cgrph.echnks(eneu)),fkern,sol1(idneu), ...
    targets(:,in2)); 
t1 = toc(start);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t1)

fkernsrc = kernel('helm','s',ks(1));
targinfo = []; targinfo.r = targets(:,in); 
uin1 = fkernsrc.fmm(1e-12,srcinfo,targinfo,charges);
utot1 = uscat1(:)-uin1(:);
utot2 = uscat2(:);

figure(3)
t = tiledlayout(1,2,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg));
zztarg(in) = uscat1;
zztarg(in2) = uscat2;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(cgrph,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg));
zztarg(in) = utot1;
zztarg(in2) = uscat2;
h=pcolor(xxtarg,yytarg,log10(abs((zztarg))));
set(h,'EdgeColor','none')
hold on
plot(cgrph,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
colorbar()

title('$\log10($error$)$','Interpreter','latex','FontSize',12)
relerr1 = max(abs(utot1));
relerr2 = max(abs(uscat2));
relerr = max(relerr1, relerr2);
fprintf('relative field error %5.2e\n',relerr);
assert(relerr < 1e-4)