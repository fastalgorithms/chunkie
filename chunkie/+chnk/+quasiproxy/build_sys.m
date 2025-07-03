function fullsys = build_sys(chnkr,kh1,kh2,d,theta,ht,hb,xleft,K)

% parameters that we may want user to have option to pass in the future
nproxy = 160;
relrad=2;
nwall=120; % number of points in each wall segment 
nwall_TD=60; % top and bottom wall for RB expansion

% bloch phase
ima = sqrt(-1);
alpha = exp(ima*kh1*d*cos(theta));

% make neighbors
chnkr_l = chnkr;
chnkr_l.r(1,:) = chnkr_l.r(1,:)-d;
chnkr_r = chnkr;
chnkr_r.r(1,:) = chnkr_r.r(1,:)+d;


%build system

[A0,AL,AR] = build_A_mat(chnkr,chnkr_l,chnkr_r,kh1,kh2);

Amat = A0 + alpha^(-1)*AL + alpha*AR;

[Bmat,Cproxy1,Cproxy2,pw] = make_B_matrix(chnkr,d,nproxy,kh1,kh2,ht,hb,relrad);
% 
% plot(Cproxy1.r(1,:),Cproxy1.r(2,:),'rx');
% plot(Cproxy2.r(1,:),Cproxy2.r(2,:),'gx')

[Cmat,lwall1,lwall2] = make_C_matrix(chnkr,kh1,kh2,alpha,ht,hb,d,nwall,xleft);

% plot(lwall1.r(1,:),lwall1.r(2,:),'k.')
% plot(lwall2.r(1,:),lwall2.r(2,:),'g.')


[Dmat] = make_D_matrix(Cproxy1,Cproxy2,pw,lwall1,lwall2,kh1,kh2,d,alpha);


[Zmat,uwall,dwall] = make_Z_matrix(chnkr,chnkr_l,chnkr_r,kh1,kh2,alpha,ht,hb,nwall_TD);

% plot(uwall.r(1,:),uwall.r(2,:),'rx');
% plot(dwall.r(1,:),dwall.r(2,:),'gx')

Vmat = make_V_matrix(Cproxy1,Cproxy2,pw,kh1,kh2,uwall,dwall);

Wmat = make_W_matrix(kh1,kh2,uwall,dwall,d,theta,K);


%tristan proposed
Cproxy = cell(3,1);
lwall = cell(2,1);
udwall = cell(2,1);


Cproxy{1} = Cproxy1;
Cproxy{2} = Cproxy2;
Cproxy{3} = pw;

lwall{1} = lwall1;
lwall{2} = lwall2;

udwall{1} = uwall;
udwall{2} = dwall;


fullsys=[];
% add matrices
fullsys.Amat = Amat;
fullsys.Bmat = Bmat;
fullsys.Cmat = Cmat;
fullsys.Dmat = Dmat;
fullsys.Vmat = Vmat;
fullsys.Zmat = Zmat;
fullsys.Wmat = Wmat;
% add proxies
fullsys.Cproxy = Cproxy;
% add geometry information
fullsys.lwall = lwall;
fullsys.udwall = udwall;

% store domain information
fullsys.d = d;
fullsys.ht = ht;
fullsys.hb = hb;

fullsys.khs = [kh1,kh2];

% build region labeler
fullsys.cgrph_lab = build_reg_labeler(chnkr,hb,ht);

return



function [A0,AL,AR] = build_A_mat(chnkr,chnkr_l,chnkr_r, kh1, kh2)
% A0 = self interaction
% AL = interaction with right neighbor
% AR = interaction with left neighbor

% define the necessary kernels
Ak1k2 = kernel('helmdiff','all',[kh1,kh2]);

% build A0
% self contributions - continuity of u and it's flux
A0 = chunkermat(chnkr,Ak1k2);
A0 = diag((-1).^(1:size(A0,2))) + A0;
% build interactions with left neighbor A:
AL = chunkerkernevalmat(chnkr_l, Ak1k2, chnkr);
% build interactions with right neighbor AR
AR = chunkerkernevalmat(chnkr_r, Ak1k2, chnkr);
return

function [Bmat,Proxy1,Proxy2,pw] = make_B_matrix(chnkr,d,nproxy,kh1,kh2,ht,hb,relrad)

% [proxy,pnorm,pw] = proxy_circ_pts(nproxy);
% 
 xc1 = [0.5*(chnkr.rstor(1,end,end)+chnkr.rstor(1,1,1));chnkr.rstor(2,1,1)];
 xc2 = [0.5*(chnkr.rstor(1,end,end)+chnkr.rstor(1,1,1));chnkr.rstor(2,1,1)];

[pr, pnorm, pw] = chnk.flam.proxy_circ_pts(nproxy); 

pw = relrad*d*pw/1.5;

Proxy1 = [];
Proxy1.r = relrad*d*pr/1.5+xc1;
Proxy1.n = pnorm;

Proxy2 = [];
Proxy2.r = relrad*d*pr/1.5+xc2;
Proxy2.n = pnorm;


B11 = construct_proxy_all(kh1, Proxy1,pw,chnkr);
B22 = construct_proxy_all(kh2, Proxy2,pw,chnkr);

Bmat =[B11 -B22];

return


function A = construct_proxy_all(kh, Proxy,pw,chnkr)

repcoef = [1,1i*kh];
A = chnk.helm2d.kern(kh,Proxy,chnkr,'c2trans',repcoef).*pw.';

return
function [Cmat,lwall1,lwall2] = make_C_matrix(chnkr,kh1,kh2,alpha,ht,hb,d,n,xval)


% make the wall for top layer
[xs,~,~,~] = lege.exps(n);

% Linear map from[-1,1] to [a,b]
yy=(chnkr.r(2,1)*(1-xs)+ht*(1+xs))/2.';

lwall1 = [];
lwall1.r = [xval*ones(1,n);yy.'];
lwall1.n = [ones(1,n);zeros(1,n)];

    rwall1=lwall1;
    rwall1.r(1,:)=rwall1.r(1,:)+d;

    [matR1,matL1] = make_Cblock(lwall1,rwall1,chnkr,kh1,d);


% make the wall for bottom
yy=(hb*(1-xs)+chnkr.r(2,1)*(1+xs))/2;
lwall2 = [];
lwall2.r = [xval*ones(1,n);yy.'];
lwall2.n = [ones(1,n);zeros(1,n)];
    rwall2=lwall2;
    rwall2.r(1,:) = rwall2.r(1,:) +d;

    [matR2,matL2] = make_Cblock(lwall2,rwall2,chnkr,kh2,d);

    Rht = [matR1; matR2];
    Lft = [matL1; matL2];

 Cmat = alpha^(-2)*Rht -alpha*Lft;

return

function [matR,matL] = make_Cblock(lwall,rwall,chnkr,kh,d)


% define the necessary kernels

lwalltmp = lwall;
lwalltmp.r(1,:,:) = lwalltmp.r(1,:,:)-d;
rwalltmp = rwall;
rwalltmp.r(1,:,:) = rwalltmp.r(1,:,:)+d;

ww = chnkr.wts; ww = ww(:).';
% duplicate weights for second density
ww = repmat(ww,2,1); ww = ww(:).';

repcoef = ones(2,2);
matR = chnk.helm2d.kern(kh,chnkr,rwalltmp,'all',repcoef).*ww;
matL = chnk.helm2d.kern(kh,chnkr,lwalltmp,'all',repcoef).*ww;

return


function [Dmat] = make_D_matrix(Proxy1,Proxy2,pw,lwall1,lwall2,kh1,kh2,d,alpha)


nwall = 2*(numel(lwall1.r)/2);
nproxy = length(pw);

rwall1 = lwall1;
rwall1.r(1,:) = rwall1.r(1,:)+d;

[DL,DR]=make_Dblock(lwall1,rwall1,Proxy1,pw,kh1);

D1 =  alpha^(-1)*DR-DL;


rwall2 = lwall2;
rwall2.r(1,:) = rwall2.r(1,:)+d;

[DL,DR]=make_Dblock(lwall2,rwall2,Proxy2,pw,kh2);
D2 =  alpha^(-1)*DR-DL;

Dmat = [D1 zeros(nwall,nproxy);...
              zeros(nwall,nproxy),D2];
return
% 
function [DL,DR]=make_Dblock(lwall,rwall,Proxy,pw,kh)
% - Evaluate operator D in Cho and Barnett (2015)

repcoef = [1,1i*kh];
DR = chnk.helm2d.kern(kh,Proxy,rwall,'c2trans',repcoef).*pw.';
DL = chnk.helm2d.kern(kh,Proxy,lwall,'c2trans',repcoef).*pw.';

return

function [Zmat,uwall,dwall] = make_Z_matrix(chnkr,chnkr_l,chnkr_r,kh1,kh2,alpha,ht,hb,nwall_TD)

x_lft=-0.5;
xtmp = (1+(2*(1:nwall_TD)/nwall_TD - 1 - 1/nwall_TD))/2;
xpts = (1-xtmp)+x_lft;

uwall = [];
uwall.r = [xpts;ht*ones(1,nwall_TD)];
uwall.n = [zeros(1,nwall_TD);ones(1,nwall_TD)];


% make top
[A0,AL,AR]=make_Zblock(uwall,chnkr,chnkr_l,chnkr_r,kh1);

Z1 = A0+alpha^(-1)*AL + alpha*AR;

% make bottom 
dwall = [];
dwall.r = [xpts;hb*ones(1,nwall_TD)];
dwall.n = [zeros(1,nwall_TD);ones(1,nwall_TD)];

[A0,AL,AR]=make_Zblock(dwall,chnkr,chnkr_l,chnkr_r,kh2);

Z2 = A0+alpha^(-1)*AL + alpha*AR;


Zmat = [Z1;Z2];
return
% 
% 
function [A0,AL,AR]=make_Zblock(xpts,C,CL,CR,kh)
% -Evaluate operator Z in Cho and Barnett (2015)

ww = C.wts; ww = ww(:).';
% duplicate weights for second density
ww = repmat(ww,2,1); ww = ww(:).';

repcoef = ones(2,2);
AL = chnk.helm2d.kern(kh,CL,xpts,'all',repcoef).*ww;
A0 = chnk.helm2d.kern(kh,C,xpts,'all',repcoef).*ww;
AR = chnk.helm2d.kern(kh,CR,xpts,'all',repcoef).*ww;

return

function Vmat = make_V_matrix(Proxy1,Proxy2,pw,kh1,kh2,uwall,dwall)

nwall_TD = size(uwall.r,2);
nproxy = length(pw);

V1 = construct_proxy_all(kh1,Proxy1,pw,uwall);
V2 = construct_proxy_all(kh2,Proxy2,pw,dwall);

Vmat = [V1 zeros(2*nwall_TD,nproxy); zeros(2*nwall_TD,nproxy) V2];
return

function Wmat = make_W_matrix(kh1,kh2,uwall,dwall,d,theta,K)

nwall_TD = size(uwall.r,2);
ima = sqrt(-1);

kappa = kh1*cos(theta)+2*pi*(-K:1:K)/d;
ktmp = ones(nwall_TD,1)*kappa;

ku = ones(nwall_TD,1)*sqrt(kh1^2-kappa.^2);
UU = uwall.r(1,:)'*ones(1,2*K+1);

W1 = zeros(2*size(UU,1),size(UU,2));
W1(1:2:end,:) = -exp(ima*UU.*ktmp);
W1(2:2:end,:) = -ima*ku.*exp(ima*UU.*ktmp);

kd = ones(nwall_TD,1)*sqrt(kh2^2-kappa.^2);
DD = dwall.r(1,:)'*ones(1,2*K+1);

W2 = zeros(2*size(DD,1),size(DD,2));
W2(1:2:end,:) = -exp(ima*DD.*ktmp);
W2(2:2:end,:) = ima*kd.*exp(ima*DD.*ktmp);

KK=2*K+1;
loc_IND1=1:2*nwall_TD;
loc_INDn=(2*nwall_TD+1):(2*nwall_TD*2);
Wmat=zeros(2*nwall_TD*2, 2*KK);

Wmat(loc_IND1,1:KK)=W1;
Wmat(loc_INDn,(KK+1):2*KK)=W2;

return


function cgrph_lab = build_reg_labeler(chnkr,hb,ht)
% detect left and right ends
vend = chunkends(chnkr); 
[~,il] = min(vend(1,:));
[~,ir] = max(vend(1,:));

lend = vend(:,il);
rend = vend(:,ir);

% build chunkgraph out of chunker, and all of the walls
verts = [lend,rend, [lend(1);ht],[rend(1);ht],[lend(1);hb],[rend(1);hb] ];
edge2vert = [[1;2],[3;4],[5;6],[1;3],[1;5],[2;4],[2;6]];

cgrph_lab = chunkgraph(verts,edge2vert);
cgrph_lab.echnks(1) = chnkr;
cgrph_lab.regions = findregions(cgrph_lab);
return