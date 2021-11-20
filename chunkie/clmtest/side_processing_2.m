targs = clmparams.xylim;
ks = clmparams.k;

kmax = max(ks);

ymin = targs(3);
ymax = targs(4);
xmax = targs(2);
xmin = targs(1);

x_used = [];
y_used = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = round((ymax-ymin)*kmax*2*pi);

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmax*ones(size(tys));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp = chnk.flam.rand_fft_transf(A12,nmax); toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmin*ones(size(tys));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

APP = A12;
A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp2 = chnk.flam.rand_fft_transf(A12,nmax); toc;

atmp = [atmp;atmp2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np = round((xmax-xmin)*kmax*2*pi);

txs = (0:np)/np*(xmax-xmin) + xmin;
tys = ymin*ones(size(txs));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp3 = chnk.flam.rand_fft_transf(A12,nmax); toc;

atmp = [atmp;atmp3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = 2*round((xmax-xmin)*kmax*2*pi);

txs = (0:np)/np*(xmax-xmin) + xmin;
tys = ymax*ones(size(txs));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp4 = chnk.flam.rand_fft_transf(A12,nmax); toc;

atmp = [atmp;atmp4];

tic; btmp = chnk.flam.rand_fft_transf(atmp,nmax); toc;

%tic; [SK,RD,T] = id(atmp,tol); toc;

atmp_k1 = atmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zk = ks(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = round((ymax-ymin)*kmax*2*pi);

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmax*ones(size(tys));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
%zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp = chnk.flam.rand_fft_transf(A12,nmax); toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmin*ones(size(tys));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
%zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

APP = A12;
A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp2 = chnk.flam.rand_fft_transf(A12,nmax); toc;

atmp = [atmp;atmp2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np = round((xmax-xmin)*kmax*2*pi);

txs = (0:np)/np*(xmax-xmin) + xmin;
tys = ymin*ones(size(txs));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
%zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp3 = chnk.flam.rand_fft_transf(A12,nmax); toc;

atmp = [atmp;atmp3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = 2*round((xmax-xmin)*kmax*2*pi);

txs = (0:np)/np*(xmax-xmin) + xmin;
tys = ymax*ones(size(txs));

x_used = [x_used,txs];
y_used = [y_used,tys];

[Xt,Xs] = meshgrid(txs,r_im(1,:));
[Yt,Ys] = meshgrid(tys,r_im(2,:));

Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
%zk = ks(1);
A1 = besselh(0,1,zk*Z);
A2 = besselh(1,1,zk*Z);

A12 = transpose([A1,A2]);
nmax = 400;
tol  = 10^(-14);
tic; atmp4 = chnk.flam.rand_fft_transf(A12,nmax); toc;

atmp = [atmp;atmp4];

tic; btmp = chnk.flam.rand_fft_transf(atmp,nmax); toc;

%tic; [SK,RD,T] = id(atmp,tol); toc;

atmp_k2 = atmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atmp_2f = [atmp_k1;atmp_k2];
tic; btmp_2f = chnk.flam.rand_fft_transf(atmp_2f,nmax); toc;
tic; [SK,RD,T] = id(btmp_2f,tol); toc;

