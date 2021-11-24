targs = clmparams.xylim;
ks = clmparams.k;
tol= 10^(-12);
kmax = max(ks);

ymin = targs(3);
ymax = targs(4);
xmax = targs(2);
xmin = targs(1);

x_used = [];
y_used = [];

tstart = tic();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = round((ymax-ymin)*kmax/(2*pi)*12);

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmax*ones(size(tys));

x_used = [x_used,txs];
y_used = [y_used,tys];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmin*ones(size(tys));

x_used = [x_used,txs];
y_used = [y_used,tys];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np = round((xmax-xmin)*kmax/(2*pi)*12);

txs = (0:np)/np*(xmax-xmin) + xmin;
tys = ymin*ones(size(txs));

x_used = [x_used,txs];
y_used = [y_used,tys];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = 2*round((xmax-xmin)*kmax/(2*pi)*12);

txs = (0:np)/np*(xmax-xmin) + xmin;
tys = ymax*ones(size(txs));

x_used = [x_used,txs];
y_used = [y_used,tys];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts_ki2r =[];
opts_ki2r.nmax  = 800;
opts_ki2r.ifcomp= true;

xs = r_im(1,:);
ys = r_im(2,:);
xt = x_used;
yt = y_used;

zk = ks(1);
[atmp] = clm.get_kern_im2re(xs,ys,xt,yt,zk,opts_ki2r);
%tic; [SK,RD,T] = id(atmp,tol); toc;
%size(SK)

zk = ks(2);
[atmp2] = clm.get_kern_im2re(xs,ys,xt,yt,zk,opts_ki2r);
%tic; [SK,RD,T] = id(atmp2,tol); toc;
%size(SK)

zk = ks(3);
[atmp3] = clm.get_kern_im2re(xs,ys,xt,yt,zk,opts_ki2r);
%tic; [SK,RD,T] = id(atmp3,tol); toc;
%size(SK)

zk = ks(4);
[atmp4] = clm.get_kern_im2re(xs,ys,xt,yt,zk,opts_ki2r);
%tic; [SK,RD,T] = id(atmp4,tol); toc;
%size(SK)

zk = ks(5);
[atmp5] = clm.get_kern_im2re(xs,ys,xt,yt,zk,opts_ki2r);
%tic; [SK,RD,T] = id(atmp5,tol); toc;
%size(SK)

atmp_5f = [atmp;atmp2;atmp3;atmp4;atmp5];
tic; btmp_5f = chnk.flam.rand_fft_transf(atmp_5f,nmax); toc;
%tic; [SK,RD,T] = id(atmp_5f,tol); toc;
size(SK)
tic; [SK,RD,T] = id(btmp_5f,tol); toc;
size(SK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp_mat = zeros([numel(SK),numel(xs)]);
exp_mat(:,SK) = eye(numel(SK));
exp_mat(:,RD) = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zk = ks(3);
x0 = rand(1)*(xmax-xmin) + xmin;
y0 = rand(1)*(ymax-ymin) + ymin;

dan = rand(size(xs));
ds0 = rand(size(xs));
ds1 = rand(size(xs));

zz = sqrt((xs-x0).^2+(ys-y0).^2);
h0 = besselh(0,1,zk*zz);
h1 = besselh(1,1,zk*zz);

h1 = h1./zz.*(cos(dan).*(xs-x0)+1i*sin(dan).*(ys-y0));

h0_eff = h0.*ds0;
h1_eff = h1.*ds1;

field_full = sum(h0_eff+h1_eff)

s0 = ds0;
s1a= ds1.*cos(dan);
s1b= 1i*ds1.*sin(dan);

s_red_0 = exp_mat*transpose(s0);
s_red_1a= exp_mat*transpose(s1a);
s_red_1b= exp_mat*transpose(s1b);

xskel = xs(SK);
yskel = ys(SK);

zskel = sqrt((xskel-x0).^2+(yskel-y0).^2);
h0_skel = besselh(0,1,zk*zskel);
h1_skel = besselh(1,1,zk*zskel);

h1_skel = h1_skel./zskel;
h1_skel_1a = h1_skel.*(xskel-x0).*transpose(s_red_1a);
h1_skel_1b = h1_skel.*(yskel-y0).*transpose(s_red_1b);

h0_skel = h0_skel.*transpose(s_red_0);

field_skel = sum(h0_skel+h1_skel_1a+h1_skel_1b)

err = field_full -field_skel

toc(tstart)
