% Defines a smooth, 1D interface and solves the corresponding 2D
% time-harmonic Klein Gordon equation on either side of the interface.
% We impose continuity and a jump in the normal derivative across the
% interface. We also impose outgoing radiation conditions at infinity 
% along the interface.


clear all
close('all')
ifsvdplot = false;

m=2;
E=1;

xmin=-60;
xmax=60;

cparams = [];
cparams.ta = xmin;
cparams.tb = xmax;
cparams.ifclosed = 0;
pref = [];
pref.nchmax=100000;

cparams.nchmin=24;

a=1.0;
b=1;
c=0.0;
d=1;
[chnkr,ab] = chunkerfunc(@(t) nonflatinterface(t,a,b,c,d),cparams,pref);

chnkr = sort(chnkr);

n = chnkr.k*chnkr.nch;

xs = chnkr.r(1,:);
xs = xs(:);
xx = xs(:);
ys = chnkr.r(2,:);
ys = ys(:);
yy = ys(:);
kh = 1j*sqrt(m^2-E^2);

x0 = 0;
y0 = 1.5;
source = [x0 y0];
rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
rr = rr(:);
u_test = (1i/4)*besselh(0,1,kh*rr);

decay_per_c = m*(xmax-xmin)/chnkr.nch;

nchpad = 2*ceil(-log(10^(-16))/abs(decay_per_c));
ich1 = nchpad;
ich2 = chnkr.nch-nchpad;
istart = (ich1-1)*chnkr.k +1;
iend = (ich2-1)*chnkr.k;

uu = -2*m*u_test(istart:iend);
[sol_gmres,vdens] = fast_solve_wrap(uu,m,E,chnkr,istart,iend,ab);

num_xts = 240;
xts_min=-6;
xts_max=6;
num_yts = 120;
yts_min=-2;
yts_max=4;
xts = linspace(xts_min,xts_max,num_xts);
yts = linspace(yts_min,yts_max,num_yts);
[XX,YY] = ndgrid(xts,yts);
XX = XX(:);
YY = YY(:);

opts = [];
opts.fac = 0.5;
flags = flagnear(chnkr,[XX.';YY.'],opts);
[rowloc,colloc] = find(flags);
rowloc = reshape(rowloc,[],1);
colloc = reshape(colloc,[],1);
zk = 1j*sqrt(m^2-E^2);
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'s',1);
targs = [XX(rowloc).';YY(rowloc).'];
if length(targs)>0
    fints = chunkerkerneval(chnkr,fkern,vdens,targs,opts);
end

wts= weights(chnkr);
wts= wts(:);

rs = chnkr.r(1:2,:);
srcinfo = [];
srcinfo.sources = rs;
srcinfo.charges = (vdens.*(wts(:))).';
ns = numel(wts);

rtargs = [XX.';YY.'];
eps = 1e-12; %1e-12
tic;
pg = 0;
pgt = 1;
zk = 1j*sqrt(m^2-E^2);
[v,~] = hfmm2d(eps,zk,srcinfo,pg,rtargs,pgt);
vpot = v.pottarg;
if length(targs)>0
    vpot(rowloc) = 0;
    vpot(rowloc) = fints;
end

XX   = reshape(XX,[numel(xts),numel(yts)]);
YY   = reshape(YY,[numel(xts),numel(yts)]);
vpot = reshape(vpot,[numel(xts),numel(yts)]);

source = [x0 y0];
rr = sqrt((XX-source(1)).^2+(YY-source(2)).^2);
u_inc = (1i/4)*besselh(0,1,kh*rr);
utot = -vpot + u_inc;

maxsolution=max(abs(vpot(:)));
figure
h=pcolor(XX,YY,real(utot)); set(h,'EdgeColor','none'); daspect([1 1 1]);
colormap(redblue)
colorbar
caxis([-maxsolution,maxsolution])





function [sol_gmres,vdens] = fast_solve_wrap(uu,m,E,chnkr,istart,iend,ab)
figs(1) = 0;

zk = 1j*sqrt(m^2-E^2);

%%%%%%%%% initialize fast apply parameters %%%%%%%%%%%%%%%

tic;

%%%% first: the single layer operator

fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'s',1);
opdims(1) = 1; opdims(2) = 1;
opts = [];
opts.nonsmoothonly = true;
tic; sys_nn = chunkermat(chnkr,fkern,opts); toc

tic; [iinds,jinds,sinds] = find(sys_nn); toc

tic;
rs = chnkr.r(:,:);
wts= weights(chnkr);
wts= wts(:);

rsrc = rs(:,jinds);
rtar = rs(:,iinds);

rdif = vecnorm(rtar-rsrc);
h0 = besselh(0,1,zk*rdif);
h1 = besselh(1,1,zk*rdif);
vcors = 0.25*1i*(h0.').*(wts(jinds));
vcors(iinds==jinds) = 0;

smat = sparse(iinds,jinds,sinds-vcors);

%%%% second: the preconditioner
%obtain parametrization of the interface.
t_coarse_1 = ab(1,:);
t_coarse_2 = ab(2,:);
chnklen=t_coarse_2 - t_coarse_1;
[xs,~] = lege.exps(chnkr.k);
tvals=zeros(chnkr.k,chnkr.nch);
for i=1:chnkr.nch
    tvals(:,i)=t_coarse_1(i)+chnklen(i)*(xs+1)/2;
end
ts=tvals(:);

%now, discretize flat interface using the above nodes.
cparamsflat = [];
cparamsflat.ta = ab(1,1);
cparamsflat.tb = ab(end,end);
cparamsflat.tsplits = [ab(1,:),ab(end,end)];
cparamsflat.ifclosed = 0;
cparamsflat.nchmin = chnkr.nch;
prefflat = [];
prefflat.nchmax = chnkr.nch;
chnkrflat = chunkerfunc(@(t) nonflatinterface(t,0,0,0,0),cparamsflat,prefflat);
chnkrflat = sort(chnkrflat);

zkE = E;
fkern2 = @(s,t) chnk.helm1d.kern(zkE,s,t,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
opts.nonsmoothonly = true;
opts.sing = 'hs'; 
start = tic; sysmat_nn = m^2*chunkermat(chnkrflat,fkern2,opts)/zkE;
t2 = toc(start)

tic;

wts= weights(chnkr);
wts= wts(:);

[iEnds,jEnds,sEnds] = find(sysmat_nn);
tEis = ts(iEnds);
tEjs = ts(jEnds);
wEts = (wts(jEnds));
vEnds = m^2*exp(1i*zkE*abs(tEis-tEjs)).*wEts/E;

sEmat = sparse(iEnds,jEnds,sEnds-vEnds);

%%%% wrap

inds = istart:iend;

fapplymat = @(u) apply_all(u,rs,wts,zk,smat,zkE,sEmat,...
    ts,m,inds);

toc

%%%%%%%%% with a fast apply


[sol_gmres,flag,relres,iter,resvec]  = gmres(fapplymat,uu,500,10^(-12),4);
if (figs(1) == 1)
    figure; plot(ts(istart:iend),log10(abs(sol_gmres)));
end

flag
relres

sol_gmres_large = zeros([numel(wts),1]);
sol_gmres_large(inds) = sol_gmres;
[vdens] = sol_gmres_large + ...
    m^2/zkE*chnk.helm1d.sweep(sol_gmres,inds,ts,wts,zkE) + ...
    sEmat*sol_gmres_large;

end





function [v] = apply_all(uin,rs,wts,zk,smat,zkE,sEmat,...
    ts,dm,inds)

vout = dm^2/zkE*chnk.helm1d.sweep(uin,inds,ts,wts,zkE);

u = zeros([numel(wts),1]);
u(inds) = uin;

charges = u + vout + sEmat*u;

%%%%%%%%%% apply single layer
srcinfo = [];
srcinfo.sources = rs;
srcinfo.charges = (charges.*(wts(:))).';
ns = numel(wts);

targ = rs;
eps = 1e-14;
tic;
pg = 1;
[v,~] = hfmm2d(eps,zk,srcinfo,pg);

v = (v.pot).'+smat*(charges);

v = charges-2*dm*v;

v = v(inds);

end



function [val] = get_gxi(xi,dm,de,dl,y0,x0,x,y)
domeg = sqrt(dm*dm-de*de);
zk = sqrt(xi^2+domeg^2);
vmat = zeros(4,4);
vmat(1,:)= [1,-exp(-zk*y0),-1,0];
vmat(2,:)= [-zk,zk*exp(-zk*y0),-zk,0];
vmat(3,:) = [0,1,exp(-zk*y0),-1];
vmat(4,:) = [0,-zk,zk*exp(-zk*y0),-zk+2*dl];

vrhs = zeros([4,1]);
vrhs(1) = 0;
vrhs(2) = 1/sqrt(2*pi);
vrhs(3) = 0;
vrhs(4) = 0;

vcoef = vmat\vrhs;
val=(double(y>y0).*vcoef(1).*exp(zk*(y0-y))...
    +double(and(y>0,y<=y0)).*(vcoef(2).*exp(-zk*y)...
    +vcoef(3).*exp(zk*(y-y0)))...
    +double(y<=0).*vcoef(4).*exp(zk*y))...
    .*exp(1i*xi*(x-x0))/sqrt(2*pi);
end