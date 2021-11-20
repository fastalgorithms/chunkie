targs = clmparams.xylim;
ks = clmparams.k;

kmax = max(ks);

ymin = targs(3);
ymax = targs(4);
xmax = targs(2);
xmin = targs(1);

np = round((ymax-ymin)*kmax*2*pi);

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmax*ones(size(tys));

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

tys = (0:np)/np*(ymax-ymin) + ymin;
txs = xmin*ones(size(tys));

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

tic; [SK,RD,T] = id(atmp,tol); toc;

AA = [APP;A12];
tic; [SK2,RD2,T2] = id(AA,tol); toc;
