function testclm4
%% This program solves the following layered medium problem.
%
%              Omega_1
%
%                 ___  ________
%                /   \/        \ 
%               /               \
%  -------------     Omega_3     ---------------------
%               \               /
%                \_____/\______/
%
%
%              Omega_2
%
%
%
%  \Delta u_i + k_i^2 u_i = 0, i=1,...,ndomain
% 
%  u_i^{in}-u_i^{ex} = f_i, i=1,...,ncurve
%
%  1/c_i^{in}du_i^{in}/dn - 1/c_i^{ex}du_i^{ex}/dn = g_i, i=1,...,ncurve
%
%  where "in" denotes interior, "ex" denotes exterior.
%
%  The unit normal vector is ALWAYS pointing outward w.r.t. the interior domain.
%
%  Features: 1. it uses complexified x-coordinates to deal with two infinite 
%  half lines; 2. it uses the RCIP method to deal with the triple junction
%  singularity; 3. the solution in each domain is represented by c D+S on
%  all four curve segments so that the results integral equations are of the
%  second kind except at two triple junctions (not a correct mathematical
%  statement, but sort of true).
% 
%
%  The representation:
%
%  u_i = sum_{j\in \Gamma} u_{i,j}, i=1,...,ndomain
%
%  u_{i,j} = c_i D_{i,\Gamma_j}[\mu_j] + S_{i,\Gamma_j}[\rho_j],
%
%  Note here that the solution is represented via a sum of layer potentials
%  on ALL curves. This is a nice trick by Leslie which reduces
%  near-hypersingular integrals to near-singular integrals.
%
%  SKIEs:
%
%  (c_i^{in}+c_i^{ex})/2 \mu_i + (c_i^{in} D_i^{ex}[\mu_i]-c_i^{ex} D_i^{in}[\mu_i])
%     + (S_i^{ex}[\rho_i]-S_i^{in}[\rho_i])
%     - \sum_{j\ne i} u_{i,j}^{in}
%     + \sum_{j\ne i} u_{i,j}^{ex}
%  = -f_i
%
%
%  (1/c_i^{in}+1/c_i^{ex})/2 \rho_i - (D'_i^{ex}[\mu_i]- D'_i^{in}[\mu_i])
%     - (1/c_i^{ex} S'_i^{ex}[\rho_i]- 1/c_i^{in} S'_i^{in}[\rho_i])
%     + 1/c_i^{in} \sum_{j\ne i} du_{i,j}^{in}/dn  
%     - 1/c_i^{ex} \sum_{j\ne i} du_{i,j}^{ex}/dn  
%  = g_i
%  
%
%
%
close all
format long e
format compact

ndomain = 3; % number of domains
ncurve = 4; % number of curve segments
chnkr(1,ncurve) = chunker();

rn = zeros(ndomain,1); 
% rn(i) is the index of refraction of the ith domain
rn(1) = 1.0;
rn(2) = 1.6;
rn(3) = 1.4;

% k0 is the wave number in vacuum
k0 = 4;

% k(i) is the wave number for the ith domain
k = k0*rn;


% coefficients in the boundary conditions on normal derivatives
coef = 1./rn.^2;

% domain indices for each curve
c = zeros(2,ncurve);
% interior domain for the ith curve
c(1,1:2) = 1; c(1,3:4) = 3;
% exterior domain for the ith curve 
c(2,1:2) = 2; c(2,3) = 1; c(2,4) = 2;


k1 = zeros(1,ncurve); % wave numbers for the interior domain
k2 = zeros(1,ncurve); % wave numbers for the exterior domain

for i=1:ncurve
  k1(i) = k(c(1,i));
  k2(i) = k(c(2,i));
end

% two circular arcs for the center eye for now
theta = zeros(1,ncurve);
% upper curve opening angle
theta(3) = pi/2.4;
% lower curve opening angle
theta(4) = pi/2.2;

% parameters for the complexification of left and right flat parts.
c1 = log(1d-2/eps)/min(k(1:2))
c2 = c1/2.5

% curve parameters
% center eye range on the real x-axis [a,b]
a = -1d1;
b = 1d1;
% length for complexification
C = 6*c2

% left flat curve [-L(1),a]
L(1) = C-a;
% right flat curve [b,L(2)]
L(2) = b+C;

% number of chunks on each curve
% should be proportional to the length*wavenumber
nch = zeros(1,ncurve);

n0 = 6;
fac = 1.2;

lambda = 2*pi/max(abs(k(1)),abs(k(2)));
nch(1) = round(fac*(a+L(1))/lambda) + n0;
nch(2) = round(fac*(L(2)-b)/lambda) + n0;

n0 = 16;
lambda = 2*pi/max(abs(k(1)),abs(k(3)));
nch(3) = round(fac*(b-a)/2/sin(theta(3)/2)/lambda) + n0;
lambda = 2*pi/max(abs(k(2)),abs(k(3)));
nch(4) = round(fac*(b-a)/2/sin(theta(4)/2)/lambda) + n0;

nch
% discretize the boundary
tab = zeros(2,ncurve);
tab(:,1) = [-L(1);a];
tab(:,2) = [b;L(2)];
for i=3:4
  tab(:,i) = [-theta(i)/2;theta(i)/2];
end

cparams = cell(1,ncurve);

for i=1:ncurve
  cparams{i}.ta = tab(1,i);
  cparams{i}.tb = tab(2,i);
  cparams{i}.ifclosed = false;
end

cpars = cell(1,ncurve);
cpars{1}.L = L(1);
cpars{1}.c1 = c1;
cpars{1}.c2 = c2;

cpars{2}.L = L(2);
cpars{2}.c1 = c1;
cpars{2}.c2 = c2;

cpars{3}.v0 = [a;0];
cpars{3}.v1 = [b;0];
cpars{3}.theta = theta(3);
cpars{3}.ifconvex = 0;

cpars{4}.v0 = [a;0];
cpars{4}.v1 = [b;0];
cpars{4}.theta = theta(4);
cpars{4}.ifconvex = 1;

% number of Gauss-Legendre nodes on each chunk
ngl = 16;
[glnodes,glwts] = lege.exps(ngl);

pref = []; 
pref.k = ngl;

% Treat the representation as if it were a 2x2 operator so that four layer 
% potentials D, S, D', S' can be evaluated together. This will avoid 
% redundant evaluation of Hankel functions.
opdims(1)=2;opdims(2)=2;

% set up GGQ machinery
[logquad] = chnk.quadggq.setuplogquad(ngl,opdims);

% log correction for self interaction in kernel-split style
isclosed = 0;
hlocal = [1];
LogC = chnk.quadjh.setuplogquad(hlocal,ngl,isclosed,opdims);
logquad.LogC  = LogC;

% define functions for curves
fcurve = cell(1,ncurve);
for icurve=1:ncurve
  fcurve{icurve} = @(t) clm.funcurve4(t,icurve,cpars{icurve});
end

% discretize the boundary
start = tic; 

for icurve=1:ncurve
  chnkr(icurve) = chunkerfuncuni(fcurve{icurve},nch(icurve),cparams{icurve},pref);
end

t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% plot geometry
figure(1)
clf
plot(chnkr,'r-','LineWidth',2)
axis equal
drawnow

% figure(2)
% clf
% quiver(chnkr)
% axis equal
% pause

rpars = [];
rpars.k = k;
rpars.c = c;
rpars.coef = coef;

isrcip = 1;

% build the system matrix
opts = [];
start = tic;
ilist=[];
[M,np,alpha1,alpha2] = clm.buildmat_fast(chnkr,rpars,opts,opdims,glwts,ilist,logquad);
if ~isrcip, M = M + eye(2*np); end
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

% compute the preconditioner R in the RCIP method for triple junctions

if isrcip
  
  start = tic;

  issymmetric = 1;
  
  opts.quad = 'jhlog';
  hlocal1 = [0.5, 0.5, 1];
  hlocal0 = [1, 0.5, 0.5];

  LogC0 = chnk.quadjh.setuplogquad(hlocal0,ngl,isclosed,opdims);
  LogC1 = chnk.quadjh.setuplogquad(hlocal1,ngl,isclosed,opdims);

  logquad.LogC0 = LogC0;
  logquad.LogC1 = LogC1;

  inds = [0, cumsum(nch)];
  % number of triple junctions
  ncorner = 2;

  nedge = 3;
  ndim = opdims(2);
  nplocal = 2*ngl*nedge*ndim;

  R = cell(1,ncorner);

  RG = speye(2*np);

  for icorner=1:ncorner
    if icorner==1
      clist = [1, 3, 4];
      isstart = [0, 0, 1];
    else
      clist = [2, 3, 4];
      isstart = [1, 1, 0];
    end
    
    if issymmetric && icorner == 2
      
    else
      rparslocal = [];
      rparslocal.k = k;
      rparslocal.c = c(:,clist);
      rparslocal.coef = coef;

      cparslocal = cell(1,nedge);
      for i=1:nedge
        cparslocal{i} = cpars{clist(i)};
        cparslocal{i}.islocal = isstart(i);
      end
      
      fcurvelocal = cell(1,nedge);
      for i=1:nedge
        fcurvelocal{i} = @(t) clm.funcurve4(t,clist(i),cparslocal{i});
      end

      [Pbc,PWbc,starL,circL,starS,circS,ilist] = rcip.setup(ngl,ndim,nedge,isstart);

      h0 = zeros(1,nedge); % chunk size at the coarsest level
      for i=1:nedge
        if isstart(i)
          h0(i) = chnkr(clist(i)).h(1);
        else
          h0(i) = chnkr(clist(i)).h(end);
        end
      end
      h0 = h0*2; % note the factor of 2 here!!!!

      nsub = 32; % level of dyadic refinement in the forward recursion for computing R

      R{icorner} = rcip.Rcomp_fast(ngl,nedge,ndim,Pbc,PWbc,nsub,...
        starL,circL,starS,circS,ilist,...
        h0,isstart,fcurvelocal,rparslocal,opts,opdims,glnodes,glwts,logquad);
    end
    
    starind = [];
    for i=1:nedge
      if isstart(i)
        starind = [starind inds(clist(i))*ngl*ndim+(1:2*ngl*ndim)];
      else
        starind = [starind inds(clist(i)+1)*ngl*ndim-fliplr(0:2*ngl*ndim-1)];
      end
    end

    M(starind,starind) = zeros(nplocal);
        
    if issymmetric && icorner == 2
      % due to pointwise block structure, need to reverse the order for
      % each edge, while keeping the same order of 2x2 blocks.
      R11 = R{1}(1:2:end,1:2:end);
      R12 = R{1}(1:2:end,2:2:end);
      R21 = R{1}(2:2:end,1:2:end);
      R22 = R{1}(2:2:end,2:2:end);
      
      indinv = [];
      n0 = 2*ngl;
      for i=1:nedge
        indinv = [indinv (i-1)*n0+(n0:-1:1)];
      end
      
      R11 = R11(indinv,indinv);
      R12 = R12(indinv,indinv);
      R21 = R21(indinv,indinv);
      R22 = R22(indinv,indinv);
      
      R2 = zeros(2*ngl*nedge*ndim);
      
      R2(1:2:end,1:2:end) = R11;
      R2(1:2:end,2:2:end) = R12;
      R2(2:2:end,1:2:end) = R21;
      R2(2:2:end,2:2:end) = R22;
      
      RG(starind,starind) = R2;
    else
      RG(starind,starind) = R{icorner};
    end
  end

  t1 = toc(start);

  fprintf('%5.2e s : time to compute R\n',t1)
end

% src(:,~i) are sources for the ith domain
nsrc = ndomain;
src = [-1, 1.3, 0; max(chnkr(3).r(2,:),[],'all')*1.3,...
  min(chnkr(4).r(2,:),[],'all')*1.3, 0];
src
hold on;plot(src(1,:),src(2,:),'r*')

% construct artificial boundary data for testing purpose
rhs = zeros(2*np,1);

for i=1:ncurve
  d1 = c(1,i); % interior domain index
  d2 = c(2,i); % exterior domain index
  
  j1 = d1 + 1; % src index for the interior domain
  if j1 > ndomain
    j1 = j1 - ndomain;
  end
  j2 = d2 + 1; % src index for the exterior domain
  if j2 > ndomain
    j2 = j2 - ndomain;
  end
  
  c1 = coef(d1);
  c2 = coef(d2);
  
  ind1 = sum(nch(1:i-1))*ngl*2+(1:2:2*nch(i)*ngl);
  ind2 = sum(nch(1:i-1))*ngl*2+(2:2:2*nch(i)*ngl);
  
  targnorm = chnk.normal2d(chnkr(i));
  nx = targnorm(1,:); nx=nx.';
  ny = targnorm(2,:); ny=ny.';
   
  [u1,gradu1]=chnk.helm2d.green(k1(i),src(:,j1),chnkr(i).r(:,:));
  du1dn = gradu1(:,1).*nx + gradu1(:,2).*ny;

  [u2,gradu2]=chnk.helm2d.green(k2(i),src(:,j2),chnkr(i).r(:,:));
  du2dn = gradu2(:,1).*nx + gradu2(:,2).*ny;
  
  rhs(ind1) = -alpha1(i)*(u1-u2);
  rhs(ind2) =  alpha2(i)*(1/c1*du1dn-1/c2*du2dn);
end

rhs = rhs(:);

% solve the linear system using gmres
start = tic; 
if isrcip
  [soltilde,it] = rcip.myGMRESR(M,RG,rhs,2*np,200,eps*20);
  sol = RG*soltilde;
  disp(['GMRES iter = ',num2str(it)])
else
  sol = gmres(M,rhs,[],1e-13,100);
end
t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% compute the solution at one target point in each domain
targ1 = [-1; 500]; % target point in domain #1
targ2 = [1.1; -100]; % target point in domain #2
targ3 = [0; 0]; % target point in domain #3

targ = [targ1, targ2, targ3]

%sol1 = sol(1:2:2*np); % double layer density
%sol2 = sol(2:2:2*np); % single layer density
%[abs(sol(1:10));abs(sol(2*ngl*(nch(1)+nch(2))+(-10:0)))]


uexact = zeros(ndomain,1);
ucomp = zeros(ndomain,1);

% compute the exact solution
for i=1:ndomain
  j=i+1;
  if j > ndomain
    j = j - ndomain;
  end
  
  uexact(i) = uexact(i) + chnk.helm2d.green(k(i),src(:,j),targ(:,i));
end

chnkrtotal = merge(chnkr);

% compute the numerical solution
for i=1:ndomain
  evalkern = @(s,t) chnk.helm2d.kern(k(i),s,t,'eval',coef(i));
  ucomp(i) = chunkerkerneval(chnkrtotal,evalkern,sol,targ(:,i));
%   skern =  @(s,t) chnk.helm2d.kern(k(i),s,t,'s',1);
%   dkern =  @(s,t) chnk.helm2d.kern(k(i),s,t,'d',1);
%   
%   for j=1:ncurve    
%     ind = sum(nch(1:j-1))*ngl+(1:nch(j)*ngl);
%     
%     dlp = chunkerkerneval(chnkr(j),dkern,sol1(ind),targ(:,i));
%     slp = chunkerkerneval(chnkr(j),skern,sol2(ind),targ(:,i));
%     
%     ucomp(i) = ucomp(i) + coef(i)*dlp + slp;
%   end
end

uerror = abs(ucomp-uexact)./abs(uexact);

%[ucomp.'; uexact.'; real(uerror)'];
fprintf('\n\n ===================\n\nerror=%d\n',norm(uerror(:)));
return

% evaluate the field in the second domain at 10000 points and record time
ngr = 30;       % field evaluation at ngr^2 points
xylim=[-10 10 -80 -50];  % computational domain
[xg,yg,targs,ngrtot] = targinit(xylim,ngr);

start = tic; 
evalkern = @(s,t) chnk.helm2d.kern(k(i),s,t,'eval',coef(2));
u = chunkerkerneval(chnkrtotal,evalkern,sol,targs);
t1 = toc(start);

fprintf('%5.2e s : time to evaluate the field at 10000 points\n',t1)

% plot out the field
fieldplot(u,1,xg,yg,xylim,ngr)

function [xg,yg,targs,ngrtot]=targinit(xylim,ngr)
xg=linspace(xylim(1),xylim(2),ngr);
yg=linspace(xylim(3),xylim(4),ngr);
ngrtot=ngr^2;
targs=zeros(2,ngrtot);
for k=1:ngr
  targs(1,(k-1)*ngr+(1:ngr)) = xg(k);
  targs(2,(k-1)*ngr+(1:ngr)) = yg; 
end

function fieldplot(u,z,xg,yg,xylim,ngr)
F1=zeros(ngr);
for k=1:ngr
  F1(1:ngr,k)=u((k-1)*ngr+(1:ngr));
end
fh =  findobj('type','figure');
nf = length(fh);
figure(nf+1)
imagesc(xg,yg,real(F1));      
colormap(jet)
axis xy
axis equal
colorbar
figure(nf+2)
imagesc(xg,yg,imag(F1));      
colormap(jet)
axis xy
axis equal
colorbar

% hold on
% np=length(z)/2;
% xy=1.1*xylim;
% zext=[xy(1);0;z(1:np);1;xy(2);xy(2)+1i*xy(3);xy(1)+1i*xy(3);xy(1)];
% fill(real(zext),imag(zext),'w','EdgeColor','w')  
% zext=[xy(2);1;z(np+1:2*np);0;xy(1);xy(1)+1i*xy(4);xy(2)+1i*xy(4);xy(2)];
% fill(real(zext),imag(zext),'w','EdgeColor','w')  
% title('Field $u({\bf x})$','Interpreter','LaTeX')
% xlabel('$x_1$','Interpreter','LaTeX')
% ylabel('$x_2$','Interpreter','LaTeX')
% axis(xylim)
% axis equal