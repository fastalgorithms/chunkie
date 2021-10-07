function testclm0
%% validate complexification of the coordinates when solving a layered 
% media problem using the chunkie code by Travis

close all
format long e
format compact

addpaths_loc();

k1=4.5; % wavenumber for domain 1
k2=4.2; % wavenumber for domain 2

% curve parameters
L = 10;
c = 8;
a = 3;

% number of chunks
nch = L*2;

% discretize domain

cparams = [];
cparams.eps = 1.0e-3;
cparams.ta = -L;
cparams.tb = L;
cparams.ifclosed = false;
cparams.nover = 0;
cparams.maxchunklen = 1.0/max(abs(k1),abs(k2)); % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16; % number of Gauss-Legendre nodes on each chunk

start = tic; 
% complexified interface
chnkr = chunkerfuncuni(@(t) clm.complexx(t,L,c,a),nch,cparams,pref); 
%chnkr = chunkerfunc(@(t) complexx(t,L,c,a),cparams,pref);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0);

% plot geometry and data
figure(1)
clf
plot(chnkr,'r.')
axis equal
%

% define kernels
sk1 =  @(s,t) chnk.helm2d.kern(k1,s,t,'s',1);
sk2 =  @(s,t) chnk.helm2d.kern(k2,s,t,'s',1);

dk1 =  @(s,t) chnk.helm2d.kern(k1,s,t,'d',1);
dk2 =  @(s,t) chnk.helm2d.kern(k2,s,t,'d',1);

spk1 = @(s,t) chnk.helm2d.kern(k1,s,t,'sprime',1);
spk2 = @(s,t) chnk.helm2d.kern(k2,s,t,'sprime',1);

dpk1 = @(s,t) chnk.helm2d.kern(k1,s,t,'dprime',1);
dpk2 = @(s,t) chnk.helm2d.kern(k2,s,t,'dprime',1);

opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
% build matrices for 8 layer potentials
S1  = 2*chunkermat(chnkr,sk1,opts);
S2  = 2*chunkermat(chnkr,sk2,opts);

D1  = 2*chunkermat(chnkr,dk1,opts);
D2  = 2*chunkermat(chnkr,dk2,opts);

S1p = 2*chunkermat(chnkr,spk1,opts);
S2p = 2*chunkermat(chnkr,spk2,opts);

D1p = 2*chunkermat(chnkr,dpk1,opts);
D2p = 2*chunkermat(chnkr,dpk2,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)
% build Kleinman–Martin equations
c = 1.22; % c is now the free parameter in Johan's paper. 
kappa = 3.0; % constant in the boundary condition on the normal derivatives

alpha1 = 1/(1+kappa);
alpha2 = kappa/(1+kappa);
alpha3 = 1/(c+kappa);
alpha4 = kappa/(c+kappa);

npts = chnkr.k*chnkr.nch; % total number of discretization points

I = eye(npts);

% system matrix, equation (25) in Johan's paper
M = [I+alpha2*D2-alpha1*D1  alpha1*(c*S2-S1);...
    -alpha4*(D2p-D1p)   I-c*alpha3*S2p+alpha4*S1p];
    
% source point for the upper half plane
src1=[0;-2];
% source point for the lower half plane
src2=[1.1;3];

% construct artificial boundary data for testing purpose
[u1,gradu1]=chnk.helm2d.green(k1,src1,chnkr.r(:,:));
[u2,gradu2]=chnk.helm2d.green(k2,src2,chnkr.r(:,:));

srcnorm = chnk.normal2d(chnkr);
nx = srcnorm(1,:); nx=nx.';
ny = srcnorm(2,:); ny=ny.';

du1dn = gradu1(:,1).*nx+gradu1(:,2).*ny;
du2dn = gradu2(:,1).*nx+gradu2(:,2).*ny;

du = -2*alpha1*(u1-u2);
dudn = 2*alpha3*(kappa*du1dn-du2dn);

rhs = [du; dudn]; rhs = rhs(:);

% solve the linear system using gmres
start = tic; sol = gmres(M,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot


% compute the solution at two target points
targ1 = [100; 1000]; % target point on the upper half plane
targ2 = [1; -1.1]; % target point on the lower half plane

sol1 = sol(1:npts); % double layer density
sol2 = sol(npts+1:end); % single layer density
figure;clf; plot(abs(sol1),'r.');hold on; plot(abs(sol2),'b-')
abs(sol1(1))
abs(sol2(1))

dt1 = chunkerkerneval(chnkr,dk1,sol1,targ1);
st1 = chunkerkerneval(chnkr,sk1,sol2,targ1);

ut1comp = dt1+st1 % u1 = D1 + S1, representation for u1 in Johan's paper
ut1exact = chnk.helm2d.green(k1,src1,targ1)
ut1error = abs(ut1comp-ut1exact)/abs(ut1exact)

dt2 = chunkerkerneval(chnkr,dk2,sol1,targ2);
st2 = chunkerkerneval(chnkr,sk2,sol2,targ2);

ut2comp = kappa*dt2 + c*st2 % representation for u2 in Johan's paper
ut2exact = chnk.helm2d.green(k2,src2,targ2)
ut2error = abs(ut2comp-ut2exact)/abs(ut2exact)
