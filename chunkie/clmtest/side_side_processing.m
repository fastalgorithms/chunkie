xvs = (1:300)/301.0*(xmax-xmin)+xmin;
yvs = (1:300)/301.0*(ymax-ymin)+ymin;
[xtt,ytt] = meshgrid(xvs,yvs);
xtt = xtt(:);
ytt = ytt(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% now find some random region

ind_r1 = find((xtt>0).*(ytt<=0));
xtt_r1 = xtt(ind_r1);
ytt_r1 = ytt(ind_r1);

[xtargs,xsrc] = meshgrid(xtt_r1,xskel);
[ytargs,ysrc] = meshgrid(ytt_r1,yskel);
xtargs = xtargs;
ytargs = ytargs;
xsrc = xsrc;
ysrc = ysrc;

zk = ks(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmax = 2*numel(xskel);

tstart = tic();

xd = xtargs - xsrc;
yd = ytargs - ysrc;
zz = sqrt(xd.^2+yd.^2);

[hh0,hh1] = hankm103(zk*zz);
%hh0 = besselh(0,1,zk*zz);
%hh1 = besselh(1,1,zk*zz);
%hh2 = besselh(2,1,zk*zz);

tic; btmp_h0 = chnk.flam.rand_fft_transf(hh0,nmax); toc;


hh1_r = hh1./zz;

%hh2_new = (2/zk)*hh1_r-hh0;
%norm(hh2-hh2_new,'fro')
hh2 = (2/zk)*hh1_r-hh0;

hh1_x = hh1_r.*xd;
hh1_y = hh1_r.*yd;

tic; btmp_h1r = chnk.flam.rand_fft_transf(hh1_r,nmax); toc;
tic; btmp_h1x = chnk.flam.rand_fft_transf(hh1_x,nmax); toc;
tic; btmp_h1y = chnk.flam.rand_fft_transf(hh1_y,nmax); toc;

hh2_r = hh2./(zz.*zz);
hh2_xx= hh2_r.*xd.*xd;
hh2_xy= hh2_r.*xd.*yd;
hh2_yy= hh2_r.*yd.*yd;

tic; btmp_h2r = chnk.flam.rand_fft_transf(hh2_r,nmax); toc;
tic; btmp_h2xx = chnk.flam.rand_fft_transf(hh2_xx,nmax); toc;
tic; btmp_h2xy = chnk.flam.rand_fft_transf(hh2_xy,nmax); toc;
tic; btmp_h2yy = chnk.flam.rand_fft_transf(hh2_yy,nmax); toc;

hh0_r = hh0./(zz.*zz);
hh0_xx= hh0_r.*xd.*xd;
hh0_xy= hh0_r.*xd.*yd;
hh0_yy= hh0_r.*yd.*yd;

tic; btmp_h0r = chnk.flam.rand_fft_transf(hh0_r,nmax); toc;
tic; btmp_h0xx = chnk.flam.rand_fft_transf(hh0_xx,nmax); toc;
tic; btmp_h0xy = chnk.flam.rand_fft_transf(hh0_xy,nmax); toc;
tic; btmp_h0yy = chnk.flam.rand_fft_transf(hh0_yy,nmax); toc;

hh1_rr= hh1./(zz.*zz.*zz);
hh1_xx= hh1_rr.*xd.*xd;
hh1_xy= hh1_rr.*xd.*yd;
hh1_yy= hh1_rr.*yd.*yd;

tic; btmp_h1rr = chnk.flam.rand_fft_transf(hh1_rr,nmax); toc;
tic; btmp_h1xx = chnk.flam.rand_fft_transf(hh1_xx,nmax); toc;
tic; btmp_h1xy = chnk.flam.rand_fft_transf(hh1_xy,nmax); toc;
tic; btmp_h1yy = chnk.flam.rand_fft_transf(hh1_yy,nmax); toc;

bmat = [btmp_h0;btmp_h1r;btmp_h1x;btmp_h1y];
size(bmat)
tic; bmat = chnk.flam.rand_fft_transf(bmat,4*nmax); toc;
bmat = [bmat;btmp_h2r;btmp_h2xx;btmp_h2xy;btmp_h2yy];
size(bmat)
tic; bmat = chnk.flam.rand_fft_transf(bmat,4*nmax); toc;
bmat = [bmat;btmp_h0r;btmp_h0xx;btmp_h0xy;btmp_h0yy];
size(bmat)
tic; bmat = chnk.flam.rand_fft_transf(bmat,4*nmax); toc;
bmat = [bmat;btmp_h1rr;btmp_h1xx;btmp_h1xy;btmp_h1yy];
size(bmat)

btmp = bmat;
%btmp2 = transpose(btmp);
tic; btmp2 = chnk.flam.rand_fft_transf(btmp,4*nmax); toc;
btmp3 = transpose(btmp2);
tic; btmp4 = chnk.flam.rand_fft_transf(btmp3,4*nmax); toc;
%btmp5 = transpose(btmp4);

%tic; btmp_red = chnk.flam.rand_fft_transf(btmp2,4*nmax); toc;
tol = 10^(-12);
tic; [SKT,RDT,~] = id(btmp4,tol); toc;



bmat_red = btmp2(SKT(:),:);
tic; [SK_fin,RD_fin,T_fin] = id(bmat_red,tol); toc;
size(SK_fin)

tdiff = toc(tstart)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% generate eval matrix

eva_mat = zeros([numel(SK_fin),numel(xtt_r1)]);
eva_mat(:,SK_fin) = eye(numel(SK_fin));
eva_mat(:,RD_fin) = T_fin;
eva_mat = transpose(eva_mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% now test

xtskel = xtt_r1(SK_fin);
ytskel = ytt_r1(SK_fin);

emax = 0;

for i=1:numel(xskel)
%xs0 = xskel(randsample(numel(xskel),1));
%ys0 = yskel(randsample(numel(yskel),1));
xs0 = xskel(i);
ys0 = yskel(i);

xdfull = xtt_r1-xs0;
ydfull = ytt_r1-ys0;
rrfull = sqrt((xdfull).^2+(ydfull).^2);

h0fu = besselh(0,1,zk*rrfull);

%%%% get reduced

xds = xtskel-xs0;
yds = ytskel-ys0;
rrs = sqrt((xds).^2+(yds).^2);

h0sk = besselh(0,1,zk*rrs);

h0_out = eva_mat*(h0sk);

ept =norm(h0_out-h0fu)/norm(h0fu);
emax = max(ept,emax);
end
disp('relative error for h0')
emax

%%%%%% now test h1 - dipole - x dir

emax = 0;

for i=1:numel(xskel)
%xs0 = xskel(randsample(numel(xskel),1));
%ys0 = yskel(randsample(numel(yskel),1));
xs0 = xskel(i);
ys0 = yskel(i);

xdfull = xtt_r1-xs0;
ydfull = ytt_r1-ys0;
rrfull = sqrt((xdfull).^2+(ydfull).^2);

h1fu = besselh(1,1,zk*rrfull);
h1fu = h1fu./rrfull.*xdfull;

%%%% get reduced

xds = xtskel-xs0;
yds = ytskel-ys0;
rrs = sqrt((xds).^2+(yds).^2);

h1sk = besselh(1,1,zk*rrs);
h1sk = h1sk./rrs.*xds;

h1_out = eva_mat*(h1sk);

ept =norm(h1_out-h1fu)/norm(h1fu);
emax = max(ept,emax);
end
disp('relative error for h1 - x')
emax

%%%%%% now test h1 - dipole - y dir
emax = 0;

for i=1:numel(xskel)
%xs0 = xskel(randsample(numel(xskel),1));
%ys0 = yskel(randsample(numel(yskel),1));
xs0 = xskel(i);
ys0 = yskel(i);

xdfull = xtt_r1-xs0;
ydfull = ytt_r1-ys0;
rrfull = sqrt((xdfull).^2+(ydfull).^2);

h1fu = besselh(1,1,zk*rrfull);
h1fu = h1fu./rrfull.*ydfull;

%%%% get reduced

xds = xtskel-xs0;
yds = ytskel-ys0;
rrs = sqrt((xds).^2+(yds).^2);

h1sk = besselh(1,1,zk*rrs);
h1sk = h1sk./rrs.*yds;

h1_out = eva_mat*(h1sk);

ept =norm(h1_out-h1fu)/norm(h1fu);
emax = max(ept,emax);
end
disp('relative error for h1 - y')
emax



%%%%%% now test h1' 

emax = 0;

for i=1:numel(xskel)
%for i=1:1
%xs0 = xskel(randsample(numel(xskel),1));
%ys0 = yskel(randsample(numel(yskel),1));
xs0 = xskel(i);
ys0 = yskel(i);

xdfull = xtt_r1-xs0;
ydfull = ytt_r1-ys0;
rrfull = sqrt((xdfull).^2+(ydfull).^2);

h0fu = besselh(0,1,zk*rrfull);
h1fu = besselh(1,1,zk*rrfull);
h2fu = besselh(2,1,zk*rrfull);

hder = zk*(h0fu-h2fu)/2./(rrfull.*rrfull) - h1fu./(rrfull.^3);
hderxx= hder.*(xdfull.^2);
hderxy= hder.*(xdfull).*(ydfull);
hderyy= hder.*(ydfull.^2);

htot = [hderxx,hderxy,hderyy];

%%%% get reduced

xds = xtskel-xs0;
yds = ytskel-ys0;
rrs = sqrt((xds).^2+(yds).^2);

h0sk = besselh(0,1,zk*rrs);
h1sk = besselh(1,1,zk*rrs);
h2sk = besselh(2,1,zk*rrs);

hd_sk = zk*(h0sk-h2sk)/2./(rrs.*rrs) - h1sk./(rrs.^3);
hdxx_sk = hd_sk.*(xds.^2);
hdxy_sk = hd_sk.*(xds).*(yds);
hdyy_sk = hd_sk.*(yds.^2);


hhs_out = eva_mat*[hdxx_sk,hdxy_sk,hdyy_sk];

ept =norm(htot-hhs_out,'fro')/norm(htot,'fro');
emax = max(ept,emax);
end
disp('relative error for hdders')
emax






