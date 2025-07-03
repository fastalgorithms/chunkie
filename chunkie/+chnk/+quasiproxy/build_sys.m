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

return



function [A0,AL,AR] = build_A_mat(chnkr,chnkr_l,chnkr_r, kh1, kh2)
% A0 = self interaction
% AL = interaction with right neighbor
% AR = interaction with left neighbor

% define the necessary kernels
Dk1 = kernel('helm', 'd', kh1);
Dk2 = kernel('helm', 'd', kh2);
Sk1 = kernel('helm', 's', kh1);
Sk2 = kernel('helm', 's', kh2);
Dpk1 = kernel('helm', 'dp', kh1);
Dpk2 = kernel('helm', 'dp', kh2);
Spk1 = kernel('helm', 'sp', kh1);
Spk2 = kernel('helm', 'sp', kh2);


% build A0
% self contributions
% continuity at interface
D = chunkermat(chnkr, Dk1)- chunkermat(chnkr, Dk2);
ntot = size(D,2);

S = chunkermat(chnkr, Sk1)- chunkermat(chnkr, Sk2);

% continuity of flux at interface
T = chunkermat(chnkr, Dpk1)- chunkermat(chnkr, Dpk2);

Sp = chunkermat(chnkr, Spk1)- chunkermat(chnkr, Spk2);

A0 = [-eye(ntot)+D S; T eye(ntot)+Sp];


% build interactions with left neighbor AL
% continuity through interface
Dl = chunkerkernevalmat(chnkr_l, Dk1, chnkr) - chunkerkernevalmat(chnkr_l, Dk2, chnkr);
Sl = chunkerkernevalmat(chnkr_l, Sk1, chnkr) - chunkerkernevalmat(chnkr_l, Sk2, chnkr);
% continuity of the flux through the interface
Tl = chunkerkernevalmat(chnkr_l, Dpk1, chnkr) - chunkerkernevalmat(chnkr_l, Dpk2, chnkr);
Spl = chunkerkernevalmat(chnkr_l, Spk1, chnkr) - chunkerkernevalmat(chnkr_l, Spk2, chnkr);
% full interaction
AL = [Dl Sl; Tl Spl];

% build interactions with left neighbor AR
% continuity through interface
Dr = chunkerkernevalmat(chnkr_r, Dk1, chnkr) - chunkerkernevalmat(chnkr_r, Dk2, chnkr);
Sr = chunkerkernevalmat(chnkr_r, Sk1, chnkr) - chunkerkernevalmat(chnkr_r, Sk2, chnkr);
% continuity of the flux through the interface
Tr = chunkerkernevalmat(chnkr_r, Dpk1, chnkr) - chunkerkernevalmat(chnkr_r, Dpk2, chnkr);
Spr = chunkerkernevalmat(chnkr_r, Spk1, chnkr) - chunkerkernevalmat(chnkr_r, Spk2, chnkr);
% full interaction
AR = [Dr Sr; Tr Spr];


clear D S T Sp Dl Sl Tl Spl Dr Sr Tr Spr
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

B11 =[chnk.quasiproxy.construct_proxy(kh1,Proxy1,pw,chnkr);...
    construct_proxy_flux(kh1,Proxy1,pw,chnkr)];
B22 = [chnk.quasiproxy.construct_proxy(kh2,Proxy2,pw,chnkr);...
    construct_proxy_flux(kh2,Proxy2,pw,chnkr)];


Bmat =[B11 -B22];


return


function A = construct_proxy_flux(kh,Proxy,pw,chnkr)

Dp = kernel('helm', 'dp', kh);
Sp = kernel('helm', 'sp', kh);

ima=sqrt(-1);

A = (Dp.eval(Proxy, chnkr)+ima*kh*Sp.eval(Proxy,chnkr)).*pw.';

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
D = kernel('helm', 'd', kh);
S = kernel('helm', 's', kh);
Dp = kernel('helm', 'dp', kh);
Sp = kernel('helm', 'sp', kh);

lwalltmp = lwall;
lwalltmp.r(1,:,:) = lwalltmp.r(1,:,:)-d;
rwalltmp = rwall;
rwalltmp.r(1,:,:) = rwalltmp.r(1,:,:)+d;

ww = chnkr.wts;
ww = reshape(ww,1,numel(ww));

DL1 = D.eval(chnkr,lwalltmp).*ww;
DR1 = D.eval(chnkr,rwalltmp).*ww;

SL1 = S.eval(chnkr,lwalltmp).*ww;
SR1 = S.eval(chnkr,rwalltmp).*ww;

TL1 = Dp.eval(chnkr,lwalltmp).*ww;
TR1 = Dp.eval(chnkr,rwalltmp).*ww;

DsL1 = Sp.eval(chnkr,lwalltmp).*ww;
DsR1 = Sp.eval(chnkr,rwalltmp).*ww;

matR=[DR1, SR1; TR1, DsR1];
matL=[DL1, SL1; TL1, DsL1];

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

function [DL,DR]=make_Dblock(lwall,rwall,Proxy,pw,kh)
% - Evaluate operator D in Cho and Barnett (2015)

% define the necessary kernels
D = kernel('helm', 'd', kh);
S = kernel('helm', 's', kh);
Dp = kernel('helm', 'dp', kh);
Sp = kernel('helm', 'sp', kh);

ima=sqrt(-1);

% (ones(nwall/2,1)*ww_proxy)
S_L = (D.eval(Proxy, lwall)+ima*kh*S.eval(Proxy,lwall)).*pw.';
S_R = (D.eval(Proxy, rwall)+ima*kh*S.eval(Proxy,rwall)).*pw.';



D_L = (Dp.eval(Proxy, lwall)+ima*kh*Sp.eval(Proxy,lwall)).*pw.';
D_R = (Dp.eval(Proxy, rwall)+ima*kh*Sp.eval(Proxy,rwall)).*pw.';

DR=[S_R; D_R];
DL=[S_L;D_L];

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


function [A0,AL,AR]=make_Zblock(xpts,C,CL,CR,kh)
% -Evaluate operator Z in Cho and Barnett (2015)

% define the necessary kernels
D = kernel('helm', 'd', kh);
S = kernel('helm', 's', kh);
Dp = kernel('helm', 'dp', kh);
Sp = kernel('helm', 'sp', kh);

ww = C.wts;
ww = reshape(ww,1,numel(ww));
wwL = CL.wts;
wwL = reshape(wwL,1,numel(wwL));
wwR = CR.wts;
wwR = reshape(wwR,1,numel(wwR));

DL = D.eval(CL,xpts).*wwL;
DD = D.eval(C,xpts).*ww;
DR = D.eval(CR,xpts).*wwR;

SL = S.eval(CL,xpts).*wwL;
SS  = S.eval(C,xpts).*ww;
SR = S.eval(CR,xpts).*wwR;

TL = Dp.eval(CL,xpts).*wwL;
T  = Dp.eval(C,xpts).*ww;
TR = Dp.eval(CR,xpts).*wwR;

DsL = Sp.eval(CL,xpts).*wwL;
Ds  = Sp.eval(C,xpts).*ww;
DsR = Sp.eval(CR,xpts).*wwR;

A0=[DD, SS;T, Ds];
AL=[DL,SL;TL, DsL];
AR=[DR,SR;TR,DsR];
return

function Vmat = make_V_matrix(Proxy1,Proxy2,pw,kh1,kh2,uwall,dwall)

nwall_TD = size(uwall.r,2);
nproxy = length(pw);

V1 = [chnk.quasiproxy.construct_proxy(kh1,Proxy1,pw,uwall);...
    construct_proxy_flux(kh1,Proxy1,pw,uwall)];
V2 = [chnk.quasiproxy.construct_proxy(kh2,Proxy2,pw,dwall);...
    construct_proxy_flux(kh2,Proxy2,pw,dwall)];


 Vmat = [V1 zeros(2*nwall_TD,nproxy); zeros(2*nwall_TD,nproxy) V2];
return

function Wmat = make_W_matrix(kh1,kh2,uwall,dwall,d,theta,K)

nwall_TD = size(uwall.r,2);
ima = sqrt(-1);

kappa = kh1*cos(theta)+2*pi*(-K:1:K)/d;
ktmp = ones(nwall_TD,1)*kappa;

ku = ones(nwall_TD,1)*sqrt(kh1^2-kappa.^2);
UU = uwall.r(1,:)'*ones(1,2*K+1);

W1 = [-exp(ima*UU.*ktmp); ...
    -ima*ku.*exp(ima*UU.*ktmp)];

kd = ones(nwall_TD,1)*sqrt(kh2^2-kappa.^2);
DD = dwall.r(1,:)'*ones(1,2*K+1);

W2 = [-exp(ima*DD.*ktmp); ...
    ima*kd.*exp(ima*DD.*ktmp)];


KK=2*K+1;
loc_IND1=1:2*nwall_TD;
loc_INDn=(2*nwall_TD+1):(2*nwall_TD*2);
Wmat=zeros(2*nwall_TD*2, 2*KK);

Wmat(loc_IND1,1:KK)=W1;
Wmat(loc_INDn,(KK+1):2*KK)=W2;

return
