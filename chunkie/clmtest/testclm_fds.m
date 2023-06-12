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
% The geometry gets prescribed through the geom_class structure. To
% see how the geometry needs to be prescribed, see the geometry_design
% doc
%
%
% We solve the linear system using gmres,
% and postprocess using a combination of an fmm and a fast direct solver
%
% This code depends on both the FLAM and fmm2d packages.
%
% For the default geometry being loaded, the computation time is about
% 30s.
%
%
close all
format long e
format compact

geom_class = clm.read_geom_clm9();
geom_class.ppw = 5;
clmparams = clm.setup_geom(geom_class);
tol = 1e-7;

chnkr = clm.get_geom_clmparams(clmparams);



% plot geometry
fontsize = 20;
figure(1)
clf
plot(chnkr,'r-','LineWidth',2)
hold on
if ~isempty(clmparams.verts)
  plot(clmparams.verts(1,:),clmparams.verts(2,:),'k.','MarkerSize',20)
end
axis equal
title('Boundary curves','Interpreter','LaTeX','FontSize',fontsize)
xlabel('$x_1$','Interpreter','LaTeX','FontSize',fontsize)
ylabel('$x_2$','Interpreter','LaTeX','FontSize',fontsize)

xlim([-10,10])
ylim([-16,4])
drawnow

src = clmparams.src_in;
targ = clmparams.src_in;

ndomain = clmparams.ndomain;
k = clmparams.k;
coef = clmparams.coef;

[targdomain,~,~] = clm.finddomain_gui(chnkr,clmparams,targ);





%% Get various operators 
start = tic; 
isrcip = 1;
[sk,~,exp_mat] = clm.get_compressed_postproc_im(chnkr,clmparams);

[Fskel1,Fskel2,skel_struct,opts_perm,M,RG] = clm.get_fds_gui(chnkr,clmparams,tol);
opts = [];

opts_rhs = [];
opts_rhs.itype = 1;
rhs = clm.get_rhs_gui(chnkr,clmparams,clmparams.npts,clmparams.alpha1, ...
    clmparams.alpha2,opts_rhs);

% solve the linear system using fds
disp(' ')
disp('Step 3: solve the linear system via GMRES.')

sol = chnk.flam.solve_2by2blk(rhs,Fskel1,Fskel2,skel_struct,opts_perm);
dt = toc(start);
disp(['Time on fds = ', num2str(dt), ' seconds'])



%% Compute exact solution
src = clmparams.src_in;
targ = clmparams.src_in;

ndomain = clmparams.ndomain;
k = clmparams.k;
coef = clmparams.coef;

[targdomain,~,~] = clm.finddomain_gui(chnkr,clmparams,targ);
[uexact,ugrad_exact] = clm.postprocess_uexact_gui(clmparams,targ, ... 
   targdomain);

chnkrtotal = merge(chnkr);


%% compute the numerical solution

list = cell(1,ndomain);
for i=1:ndomain
    list{i} = find(targdomain==i);
end

% Compute skelentonization, expansion, and eval matrices

opts_postproc_im = [];
[eva_mats,sktarg] = clm.get_evamat_postproc_im(chnkr,clmparams,targ, ...
      targdomain,sk,tol,opts_postproc_im);
[ucomp,ugrad] = clm.postprocess_sol_gui_fmm_fds(chnkr,clmparams, ...
    targ,targdomain,tol,sol,sk,exp_mat,eva_mats,sktarg);

ucomp = ucomp.';
uerror = abs(ucomp(:)-uexact(:))./abs(uexact);
disp(' ')
disp('Now check the accuracy of numerical solutions')
disp('Exact value               Numerical value           Error')  
fprintf('%0.15e     %0.15e     %7.1e\n', [real(uexact).'; real(ucomp).'; real(uerror)'])

return



uerror2 = abs(ugrad-ugrad_exact)./abs(ugrad_exact);
disp(' ')
disp('Now check the accuracy of gradient numerical solutions')
disp('Exact value               Numerical value           Error')  

aa = real(ugrad_exact(2,:));
bb = real(ugrad(2,:));
cc = uerror2(2,:);

fprintf('%0.15e     %0.15e     %7.1e\n', [aa; bb; cc])

%% Evaluate solution at a grid of targets


% evaluate the field in the second domain at 10000 points and record time
ngr = clmparams.ngr;% field evaluation at ngr^2 points
xylim=clmparams.xylim;  % computational domain
[xg,yg,targs,ntarg] = clm.targinit(xylim,ngr);

disp(' ')
disp(['Evaluate the field at ', num2str(ntarg), ' points'])
disp('Step 1: identify the domain for each point')
tic, [targdomain,tid] = clm.finddomain_gui(chnkr,clmparams,targs); toc;
tid = unique(tid);
ntid = setdiff(1:ntarg,tid);
list = cell(1,ndomain);
for i=1:ndomain
    list{i} = find(targdomain==i);
end

disp('Step 2: evaluate the total field for point sources')
[eva_mats,sktarg] = clm.get_evamat_postproc_im(chnkr,clmparams,targs, ...
      targdomain,sk,tol,opts_postproc_im);
start=tic; [ucomp,ugrad] = clm.postprocess_sol_gui_fmm_fds(chnkr,clmparams, ...
    targs,targdomain,tol,sol,sk,exp_mat,eva_mats,sktarg); dt = toc(start);
disp(['Evaluation time = ', num2str(dt), ' seconds'])

uexact = clm.postprocess_uexact_gui(clmparams,targs,targdomain);

% plot out the field
clm.fieldplot(ucomp,chnkr,xg,yg,xylim,ngr,fontsize)
title('Numerical solution','Interpreter','LaTeX','FontSize',fontsize)
clm.fieldplot(uexact,chnkr,xg,yg,xylim,ngr,fontsize)
title('Exact solution','Interpreter','LaTeX','FontSize',fontsize);


%% Solve with plane wave data

alpha = pi/4;
opts_rhs = [];
opts_rhs.itype = 2;
opts_rhs.alpha = alpha;



disp(' ')
disp(['Now calculate the field when the incident wave is a plane wave'])
disp(['incident angle = ', num2str(opts_rhs.alpha)])
rhs = clm.get_rhs_gui(chnkr,clmparams,clmparams.npts,clmparams.alpha1, ...
   clmparams.alpha2,opts_rhs);

% solve the linear system using fds
sol = chnk.flam.solve_2by2blk(rhs,Fskel1,Fskel2,skel_struct,opts_perm);
disp('Step 4: evaluate the total field for incident plane wve')
start = tic;
[u,gradu] = clm.postprocess_sol_gui_fmm_fds(chnkr,clmparams,targs, ...
       targdomain,tol,sol,sk,exp_mat,eva_mats,sktarg);
   
[~,m] = size(targs);
uinc = zeros(m,1);
uincgrad = zeros(2,m);
targlist = cell(1,ndomain);
for i=1:ndomain
    targlist{i} = find(targdomain==i);
end

idomup = find(clmparams.is_inf == 1);
idomdown = find(clmparams.is_inf == -1);

for i=1:ndomain
  if ~isempty(list{i})
    if i==idomup || i==idomdown
        [uinc(targlist{i}),gtmp] = clm.planewavetotal_gui(k(idomup), ...
           alpha,k(idomdown),targs(:,list{i}),clmparams.is_inf(i), ...
           idomup,idomdown,coef);
        uincgrad(:,targlist{i}) = gtmp.';
    end
  end
end
uinc = reshape(uinc,[1,m]);
u = u + uinc;
gradu = gradu + uincgrad;
dt = toc(start);


disp(['Evaluation time = ', num2str(dt), ' seconds'])
%fprintf('%5.2e s : time to evaluate the field at 10000 points\n',t1)

% plot out the field
clm.fieldplot(u,chnkr,xg,yg,xylim,ngr,fontsize)
title('Total field, incident plane wave angle = $\frac{3\pi}{4}$', ... 
   'Interpreter','LaTeX','FontSize',fontsize)

