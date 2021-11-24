
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
clear
format long e
format compact

% obtain physical and geometric parameters
icase = 6;
opts = [];
clmparams = clm.setup(icase,opts);

if isfield(clmparams,'k')
  k = clmparams.k;
end
if isfield(clmparams,'coef')
  coef = clmparams.coef;
end
if isfield(clmparams,'cpars')
  cpars = clmparams.cpars;
end
if isfield(clmparams,'cparams')
  cparams = clmparams.cparams;
end
if isfield(clmparams,'ncurve')
  ncurve = clmparams.ncurve;
end
if isfield(clmparams,'ndomain')
  ndomain = clmparams.ndomain;
end

vert = [];
if isfield(clmparams,'vert')
  vert = clmparams.vert;
end

if isfield(clmparams, 'nch')
  nch = clmparams.nch;
end

if isfield(clmparams, 'src')
  src = clmparams.src;
end

chnkr(1,ncurve) = chunker();

% define functions for curves
fcurve = cell(1,ncurve);
for icurve=1:ncurve
  fcurve{icurve} = @(t) clm.funcurve(t,icurve,cpars{icurve},icase);
end

% number of Gauss-Legendre nodes on each chunk
ngl = 16;

pref = []; 
pref.k = ngl;
disp(['Total number of unknowns = ',num2str(sum(nch)*ngl*2)])

% discretize the boundary
start = tic; 

for icurve=1:ncurve
  chnkr(icurve) = chunkerfuncuni(fcurve{icurve},nch(icurve), ...
    cparams{icurve},pref);
end


t1 = toc(start);

[chnkr,clmparams] = clm.get_geom_gui(icase,opts);
chnkrtotal = merge(chnkr);

%fprintf('%5.2e s : time to build geo\n',t1)

if (1 == 0)
% plot geometry
fontsize = 20;
figure(1)
clf
plot(chnkr,'r-','LineWidth',2)
hold on
if ~isempty(vert)
  plot(vert(1,:),vert(2,:),'k.','MarkerSize',20)
end
axis equal
title('Boundary curves','Interpreter','LaTeX','FontSize',fontsize)
xlabel('$x_1$','Interpreter','LaTeX','FontSize',fontsize)
ylabel('$x_2$','Interpreter','LaTeX','FontSize',fontsize)
drawnow
x = clm.get_region_pts_gui(chnkr,clmparams,3);
fill(x(1,:),x(2,:),'g')
end

% figure(2)
% clf
% quiver(chnkr)
% axis equal

isrcip = 1;
opts = [];
opts.nonsmoothonly = true;
[np,alpha1,alpha2] = clm.get_alphas_np_gui(chnkr,clmparams);
[M,RG] = clm.get_mat_gui_clm_cases(chnkr,clmparams,icase,opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 2*np;

rs = chnkrtotal.r(:,:);
rs = repelem(rs,1,2);
i_real = find((abs(imag(rs(1,:)))+abs(imag(rs(2,:)))) == 0);
i_imag = find((abs(imag(rs(1,:)))+abs(imag(rs(2,:)))) > 0);
iperm = [i_real,i_imag];

opts_perm = [];
%iperm = (2*np):(-1):1;
%iperm = randperm(nn);
%iperm = 1:nn;
opts_perm.iperm = iperm;
invperm = 1:nn;
invperm(iperm) = 1:nn;
opts_perm.invperm = invperm;
opts_perm.ns = [numel(i_real),numel(i_imag)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



opts_rhs = [];
opts_rhs.itype = 1;
rhs = clm.get_rhs_gui_clm_cases(chnkr,clmparams,np,alpha1,alpha2,opts_rhs);

if (1==0)
tprox = 2*pi*(0:100)/101;
rprox = 1.5;
prxy_pts = rprox*[cos(tprox);sin(tprox)];
rprox2= 1.4;
prxy_pts2= rprox2*[cos(tprox);sin(tprox)];
prxy_pts = [prxy_pts,prxy_pts2];

prxy_ks  = clmparams.k;

end

opdims(1) = 2;
opdims(2) = 2;
disp(np)
nn = 2*np;
M = M + speye(nn);
allt1 = @(s,t) chnk.helm2d.kern(1.0,s,t,'trans1',1);
wts = weights(chnkrtotal);
matfun = @(i,j) chnk.flam.kernbyindex_ell2(i,j,chnkrtotal,wts,allt1,...
    opdims,M,opts_perm);

% npxy = 100;
% [pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy);
% ifaddtrans = true;
% pxyfun = @(x,slf,nbr,l,ctr) chnk.flam.proxyfun(slf,nbr,l,ctr,chnkr,wts, ...
%        kern,opdims,pr,ptau,pw,pin,ifaddtrans);
xflam = chnkrtotal.r(:,:);
xflam = repelem(xflam,1,2);
xflam = xflam(:,opts_perm.iperm);

 rank_or_tol = 1e-9;

 pxyfun = [];
  opts_perm.n_offset = opts_perm.ns(1);
 matfun2 = @(i,j) chnk.flam.kernbyindex_ell2(i,j,chnkrtotal,wts,allt1,opdims,M,...
 opts_perm);
 irange = (opts_perm.ns(1)+1):(opts_perm.ns(1)+opts_perm.ns(2));
 xflam_i = [real(xflam(1,irange));imag(xflam(1,irange))];
 %xflam_i = xflam(:,irange);
 occ = 400;
 tic, Fskel2 = rskelf(matfun2,xflam_i,occ,rank_or_tol,pxyfun,opts); toc;
 
 rank_or_tol = 1e-9;
 occ = 400;
 opts = [];
 opts.lvlmax = 10;
 opts.verb = 1;
% pxyfun_rere = @(x,slf,nbr,l,ctr) pxy_fun_multi(x,slf,nbr,l,ctr,prxy_pts,prxy_ks);
 tic, Fskel = rskelf(matfun,xflam(:,1:opts_perm.ns(1)),occ,rank_or_tol,pxyfun,opts); toc;


 
%tic, M1 = matfun(1:2*np,1:2*np); toc;
opts_perm.n_offset = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tic, M1 = matfun(1:nn,1:nn); toc;
n1 = opts_perm.ns(1);
n2 = opts_perm.ns(2);

A12 = matfun(1:n1,(n1+1):nn);
A21 = matfun((n1+1):nn,1:n1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmax = 400;
tol  = 10^(-12);
[skel_struct] = chnk.flam.skel_2by2blk(Fskel,Fskel2,A12,A21,nmax,tol);


%M1 = matfun(1:nn,1:nn); toc;

if (1 == 0)
% 
% solve the linear system using gmres
disp(' ')
disp('Step 3: solve the linear system via GMRES.')
start = tic; 
if isrcip
  [soltilde,it] = rcip.myGMRESR(M,RG,rhs,2*np,1000,eps*20);
  sol = RG*soltilde;
  disp(['GMRES iterations = ',num2str(it)])
else
  sol = gmres(M,rhs,[],1e-13,1000);
end
dt = toc(start);
disp(['Time on GMRES = ', num2str(dt), ' seconds'])
end

%disp('solving')
%size(M1)
wts_lng = transpose([wts(:),wts(:)]);
size(wts_lng)
wts_lng = wts_lng(:);

%M2 = bsxfun(@times,M1,(1./sqrt(wts_lng)).');
%M2 = bsxfun(@times,M2,sqrt(wts_lng));

%matfun = @(i,j) chnk.flam.kernbyindex_ell2(i,j,chnkrtotal,wts,allt1,...
%    opdims,M,opts_perm);
%M3 = matfun(1:nn,1:nn); toc;

rhs_rescale = rhs.*sqrt(wts_lng);

%tic; sol = M3\rhs_rescale(iperm); toc;


%sol = rskelf_sv(F,rhs(iperm));
[sol] = chnk.flam.solve_2by2blk(rhs_rescale,Fskel,Fskel2,skel_struct,opts_perm);
%sol = sol(invperm)./sqrt(wts_lng);
sol = sol./sqrt(wts_lng);

%fprintf('%5.2f seconds : time for dense gmres\n',t1)

% compute the solution at one target point in each domain
% targ1 = [-1; 500]; % target point in domain #1
% targ2 = [1.1; -100]; % target point in domain #2
% targ3 = [0; 0]; % target point in domain #3
% targ = [targ1, targ2, targ3]
targ = src;
%sol1 = sol(1:2:2*np); % double layer density
%sol2 = sol(2:2:2*np); % single layer density
%abs(sol(1))


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



% compute the numerical solution
for i=1:ndomain
  evalkern = @(s,t) chnk.helm2d.kern(k(i),s,t,'eval',coef(i));
  ucomp(i) = chunkerkerneval(chnkrtotal,evalkern,sol,targ(:,i));
end

uerror = abs(ucomp-uexact)./abs(uexact);
disp(' ')
disp('Now check the accuracy of numerical solutions')
disp('Exact value               Numerical value           Error')  
fprintf('%0.15e     %0.15e     %7.1e\n', [real(uexact).'; real(ucomp).'; real(uerror)'])
return

% evaluate the field in the second domain at 10000 points and record time
ngr = 220;       % field evaluation at ngr^2 points
xylim=[-8 8 -12 4];  % computational domain
[xg,yg,targs,ntarg] = clm.targinit(xylim,ngr);

disp(' ')
disp(['Evaluate the field at ', num2str(ntarg), ' points'])
disp('Step 1: identify the domain for each point')
clist = clmparams.clist;
targdomain = clm.finddomain(chnkr,clist,targs,icase);
list = cell(1,ndomain);
for i=1:ndomain
    list{i} = find(targdomain==i);
end

if 1==1
    disp('Step 2: evaluate the total field for point sources')

    eps0 = 1e-7;
    start = tic;
    uexact = clm.postprocess_uexact_gui(clmparams,targs,targdomain);
    u = clm.postprocess_sol_gui(chnkr,clmparams,targs,targdomain,eps0,sol);
    dt = toc(start);


    disp(['Evaluation time = ', num2str(dt), ' seconds'])
    %fprintf('%5.2e s : time to evaluate the field at 10000 points\n',t1)

    % plot out the field
    clm.fieldplot(u,chnkr,xg,yg,xylim,ngr,fontsize)
    title('Numerical solution','Interpreter','LaTeX','FontSize',fontsize)
    clm.fieldplot(uexact,chnkr,xg,yg,xylim,ngr,fontsize)
    title('Exact solution','Interpreter','LaTeX','FontSize',fontsize)
end

% the incident wave is a plane wave

alpha = 3*pi/4;
opts_rhs = [];
opts_rhs.itype = 2;
opts_rhs.alpha = alpha;



disp(' ')
disp(['Now calculate the field when the incident wave is a plane wave'])
disp(['incident angle = ', num2str(opts_rhs.alpha)])
rhs = clm.get_rhs_gui(chnkr,clmparams,np,alpha1,alpha2,opts_rhs);

% solve the linear system using gmres
disp(' ')
disp('Step 3: solve the linear system via GMRES.')
start = tic; 
if isrcip
  [soltilde,it] = rcip.myGMRESR(M,RG,rhs,2*np,1000,eps*20);
  sol = RG*soltilde;
  disp(['GMRES iterations = ',num2str(it)])
else
  sol = gmres(M,rhs,[],1e-13,100);
end
dt = toc(start);
disp(['Time on GMRES = ', num2str(dt), ' seconds'])

disp('Step 4: evaluate the total field for incident plane wve')
eps0 = 1e-7;
start = tic;
u = clm.postprocess_sol_gui(chnkr,clmparams,targs,targdomain,eps0,sol);

start = tic; 
for i=1:ndomain
  if ~isempty(list{i})
    if i==1 || i==2
        u(list{i}) = u(list{i}) + clm.planewavetotal(k(1),alpha,k(2),targs(:,list{i}),i,coef);
    end
  end
end

dt = toc(start);


disp(['Evaluation time = ', num2str(dt), ' seconds'])
%fprintf('%5.2e s : time to evaluate the field at 10000 points\n',t1)

% plot out the field
clm.fieldplot(u,chnkr,xg,yg,xylim,ngr,fontsize)
title('Total field, incident plane wave angle = $\frac{3\pi}{4}$','Interpreter','LaTeX','FontSize',fontsize)


