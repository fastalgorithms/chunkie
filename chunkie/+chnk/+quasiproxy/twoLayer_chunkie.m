function twoLayer_chunkie

% addpath(genpath('/Users/gillmana/Documents/Research/Flatiron/chunkie/'))
% addpath(genpath('/Users/gillmana/Documents/Research/Flatiron/chunkie/+chnk/+flam'))

% Set wave speed in each layer
kh1 = 10;
kh2 = 5;
% angle 
theta = -pi/5;

% test target location for convergence.
xxtrg = [];
xxtrg.r = [0.2; 0.5];
xxtrg.n = [0; 0];


% parameters for the periodizing method
Npan = 40;
nproxy = 160;
relrad=2;
nwall=120; % number of points in each wall segment 
nwall_TD=60; % top and bottom wall for RB expansion
d = 1; % length of unit cell

% top and bottom
ht = 1;
hb = -1;

% look for woods anomaly
K = 20;
% number of terms in all Bragg expansions
KK=2*K+1;
kappa = kh1*cos(theta)+2*pi*[-K:1:K]/d;
ku = sqrt(kh1^2-kappa.^2);
kd = sqrt(kh2^2-kappa.^2);

tol_wood=1e-14;
if min(abs(ku))<tol_wood
    fprintf('\n Warning: Wood anormaly is detected at top layer! (checking tolerance set to %5.2e)',...
        tol_wood);
end
if min(abs(kd))<tol_wood
    fprintf('\n Warning: Wood anormaly is detected at bottom layer! (checking tolerance set to %5.2e)',...
        tol_wood);
end

ima = sqrt(-1);
alpha = exp(ima*kh1*d*cos(theta));


% Initialize the geometry
cparams = [];
cparams.ta = -0.5;
cparams.tb = 0.5;
chnkr = chunkerfuncuni(@(t) fcurve(t), Npan, cparams);

cparams_r = [];
cparams_r.ta = 0.5;
cparams_r.tb = 1.5;
chnkr_r = chunkerfuncuni(@(t) fcurve(t), Npan, cparams_r);

cparams_l = [];
cparams_l.ta = -1.5;
cparams_l.tb = -0.5;
chnkr_l = chunkerfuncuni(@(t) fcurve(t), Npan, cparams_l);

plot(chnkr, 'k.'); axis equal; hold on;
plot(chnkr_r, 'b.');
plot(chnkr_l, 'r.');

[A0,AL,AR] = build_A_mat(chnkr,chnkr_l,chnkr_r,kh1,kh2);

A = A0 + alpha^(-1)*AL + alpha*AR;

[Bmat,Cproxy1,Cproxy2,pw] = make_B_matrix(chnkr,d,nproxy,kh1,kh2,ht,hb,relrad);
% 
plot(Cproxy1.r(1,:),Cproxy1.r(2,:),'rx');
plot(Cproxy2.r(1,:),Cproxy2.r(2,:),'gx')

[Cmat,lwall1,lwall2] = make_C_matrix(chnkr,kh1,kh2,alpha,ht,hb,d,nwall,cparams.ta);

plot(lwall1.r(1,:),lwall1.r(2,:),'k.')
plot(lwall2.r(1,:),lwall2.r(2,:),'g.')


[Dmat] = make_D_matrix(Cproxy1,Cproxy2,pw,lwall1,lwall2,kh1,kh2,d,alpha);


[Zmat,uwall,dwall] = make_Z_matrix(chnkr,chnkr_l,chnkr_r,kh1,kh2,alpha,ht,hb,nwall_TD);

plot(uwall.r(1,:),uwall.r(2,:),'rx');
plot(dwall.r(1,:),dwall.r(2,:),'gx')

Vmat = make_V_matrix(Cproxy1,Cproxy2,pw,kh1,kh2,uwall,dwall);

Wmat = make_W_matrix(kh1,kh2,uwall,dwall,d,theta,K);

% make RHS (onlly nonzero part)
rhs = make_rhs(chnkr,kh1,theta);

% Blocksolve the system.
invA = inv(A);
% 
AinvB = invA*Bmat;



Schur_top=Dmat-Cmat*AinvB;
Schur_bot=Vmat-Zmat*AinvB;

num_layer = 2;
Schur3=[Schur_top, zeros(num_layer*nwall*2,2*KK);...
    Schur_bot, Wmat];

[SU,SS,SV] = svd(Schur3);

ind = find(diag(SS)>1e-13);
dia = diag(SS);
dia = dia(ind);

SU = SU(:,ind);
SS = diag(1./dia);
SV = SV(:,ind);

% A^-1*f
Ainvf = invA*rhs;

rhs_schur3=-[Cmat*Ainvf;...
    Zmat*Ainvf];

% proxy circle charges
coeff=SV*(SS*(SU'*(rhs_schur3)));

% Rayleigh Bragg coefficients (top and bottom)
bragg_a=coeff(num_layer*nproxy+1:num_layer*nproxy+2*KK);
% proxy circle charges
prox_dens=coeff(1:num_layer*nproxy);
% Interface densities
interface_dens=Ainvf-AinvB*prox_dens;

%% flux error check
    fprintf('\n All wavespeed are real, energy is preserved: ');
    KK=2*K+1;
    
    k1x=kh1*cos(theta);
    k1y=kh1*sin(theta);
        [flux_up,flux_down]=compute_flux(bragg_a(1:KK),bragg_a(KK+1:2*KK),...
        k1x,kh1,kh2,d,K);
    total_flux_up        = sum(real(flux_up));
    total_flux_down      = sum(real(flux_down));
    flux_error           = abs((total_flux_up+total_flux_down-abs(k1y))/abs(k1y));
    fprintf('Flux error est.    = %g\n\n',      flux_error);
    fprintf('\n  theta value , flux error est.    = %g',   flux_error);



% extract the different parts of the solution to the linear system.
ntot = size(A,1)/2;
 sig=interface_dens(1:ntot);
 tau=interface_dens(ntot+1:2*ntot);
 c1=prox_dens(1:nproxy);


 %  create the approximate solution at the target point
D = kernel('helm', 'd', kh1);
S = kernel('helm', 's', kh1);


ww = chnkr.wts;
ww = reshape(ww,1,numel(ww));
wwL = chnkr_l.wts;
wwL = reshape(wwL,1,numel(wwL));
wwR = chnkr_r.wts;
wwR = reshape(wwR,1,numel(wwR));

DL = D.eval(chnkr_l,xxtrg).*wwL;
DD = D.eval(chnkr,xxtrg).*ww;
DR = D.eval(chnkr_r,xxtrg).*wwR;

SL = S.eval(chnkr_l,xxtrg).*wwL;
SS  = S.eval(chnkr,xxtrg).*ww;
SR = S.eval(chnkr_r,xxtrg).*wwR;

us = alpha^(-1)*SL*tau+alpha*SR*tau+SS*tau;

ud = alpha^(-1)*DL*sig+alpha*DR*sig+DD*sig;


uproxy = construct_proxy(kh1,Cproxy1,pw,xxtrg)*c1;


uapp = us+ud+uproxy;

uapp


% fprintf('\n Npan=%d, Ngau=%d, ',...
%     Npan, Ngau);
fprintf('width of a period =%4.1e.\n', d);
fprintf('\n wave number of top layer =%5.2f.\n', kh1);
fprintf('\n Solving for single incident angle (top layer), theta=%5.2f*pi.\n',theta/pi);
fprintf('parameters: top and bottom wall=%d points,\n', nwall_TD);
fprintf('            left and right wall=%d points,\n', nwall/2);

fprintf('            proxy points per layer=%d with rel radius =%8.5e,\n',...
    nproxy, relrad);
fprintf('            num of term in expansion =2K+1, with K=%d.\n',K);
fprintf('            u_test = %d.\n',uapp);



return

% geometry
function r = fcurve(t)
    t1 = t(:);
    r = zeros(2,length(t1));
    r(1,:) = t1.';
    r(2,:) = 0.25*sin(2*pi*(t1+0.5));
    r = reshape(r, [2, size(t)]);
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

B11 =[construct_proxy(kh1,Proxy1,pw,chnkr);...
    construct_proxy_flux(kh1,Proxy1,pw,chnkr)];
B22 = [construct_proxy(kh2,Proxy2,pw,chnkr);...
    construct_proxy_flux(kh2,Proxy2,pw,chnkr)];


Bmat =[B11 -B22];


return

function A = construct_proxy(kh, Proxy,pw,chnkr)

D = kernel('helm', 'd', kh);
S = kernel('helm', 's', kh);

ima=sqrt(-1);

A = (D.eval(Proxy, chnkr)+ima*kh*S.eval(Proxy,chnkr)).*pw.';

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

% function r = vertwall_function(t,xcoor,yb,yt)
% % t parameter
% % xcoor = xcoordinate of the wall
% % yb = bottom y coordinate
% % yt = top y coordinate
% 
%     m = (yb-yt);
%     t1 = t(:);
%     r = zeros(2,length(t1));
%     r(1,:) = xcoor*ones(size(t1)).';
%     r(2,:) = m*t1+yt;
%     r = reshape(r, [2, size(t)]);
% return

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

V1 = [construct_proxy(kh1,Proxy1,pw,uwall);...
    construct_proxy_flux(kh1,Proxy1,pw,uwall)];
V2 = [construct_proxy(kh2,Proxy2,pw,dwall);...
    construct_proxy_flux(kh2,Proxy2,pw,dwall)];


 Vmat = [V1 zeros(2*nwall_TD,nproxy); zeros(2*nwall_TD,nproxy) V2];
return

function Wmat = make_W_matrix(kh1,kh2,uwall,dwall,d,theta,K)

nwall_TD = size(uwall.r,2);
ima = sqrt(-1);

kappa = kh1*cos(theta)+2*pi*[-K:1:K]/d;
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

function rhs = make_rhs(chnkr,kh1,theta)

    [u_inc,un_inc] = make_incident(theta,kh1,chnkr);
    rhs=[-u_inc;-un_inc];

return

function [u_inc,un_inc] = make_incident(theta,kh,C)
% - Incident wave generator (for the top interface)
% INPUT: 
%       theta, incident angle
%       kh, wave number of the top layer
%       C, the parameterized geometry for the top interface

ima = sqrt(-1);
nn1 = C.n(1,:).';
nn2 = C.n(2,:).';


 u_inc = exp(ima*kh*(cos(theta)*C.r(1,:)'+sin(theta)*C.r(2,:)'));
 
 un_inc = ima*kh*(nn1*cos(theta)+nn2*sin(theta)).*u_inc;


return

function [flux_up, flux_down] = compute_flux(au, ad, k1x, k1, kN, d, BR)
% - Flux error estimator, from Cho and Barnett's solver


% flux_up = 0;
% flux_down = 0;

flux_up = zeros(2*BR+1,1);
flux_down = zeros(2*BR+1,1);
for J=-BR:BR
    kappa_n  = k1x+2*pi*J/(d);
    ku_n     = sqrt(k1^2-kappa_n^2);
    kd_n     = sqrt(kN^2-kappa_n^2);
    
    flux_up(J+BR+1) = ku_n*(abs(au(J+BR+1))^2);
    flux_down(J+BR+1)=kd_n*(abs(ad(J+BR+1))^2);
%     flux_up  =  flux_up+ku_n*(abs(au(J+BR+1))^2);
%     flux_down = flux_down+kd_n*(abs(ad(J+BR+1))^2);
end

