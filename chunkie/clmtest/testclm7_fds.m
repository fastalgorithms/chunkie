
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


geom_class = clm.read_geom_clm8();
clmparams = clm.setup_geom(geom_class);


chnkr = clm.get_geom_clmparams(clmparams);

%fprintf('%5.2e s : time to build geo\n',t1)

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

x = clm.get_region_pts_gui(chnkr,clmparams,2);
plot(x(1,:),x(2,:),'g-','LineWidth',3)
% if(~isempty(xdom))
%     plot(xdom(1,:),xdom(2,:),'b-','LineWidth',3);
% end
xlim(clmparams.xylim(1:2))
ylim(clmparams.xylim(3:4))
drawnow
return


% figure(2)
% clf
% quiver(chnkr)
% axis equal

chnkrtotal = merge(chnkr);
[np,alpha1,alpha2] = clm.get_alphas_np_gui(chnkr,clmparams);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts_rhs = [];
opts_rhs.itype = 1;
rhs = clm.get_rhs_gui(chnkr,clmparams,clmparams.npts,clmparams.alpha1,clmparams.alpha2,opts_rhs);

disp(np)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%disp('solving')
%size(M1)
%tic; sol = M1\rhs(iperm); toc;
eps = 0.5e-5;
[Fskel,Fskel2,skel_struct,opts_perm,M,RG] = clm.get_fds_gui(chnkr,clmparams,eps);


%sol = rskelf_sv(F,rhs(iperm));
[sol] = chnk.flam.solve_2by2blk(rhs,Fskel,Fskel2,skel_struct,opts_perm);
%sol = sol(invperm);

%fprintf('%5.2f seconds : time for dense gmres\n',t1)

% compute the solution at one target point in each domain
% targ1 = [-1; 500]; % target point in domain #1
% targ2 = [1.1; -100]; % target point in domain #2
% targ3 = [0; 0]; % target point in domain #3
% targ = [targ1, targ2, targ3]
targ = clmparams.src_in;
%sol1 = sol(1:2:2*np); % double layer density
%sol2 = sol(2:2:2*np); % single layer density
%abs(sol(1))

ndomain = clmparams.ndomain;
uexact = zeros(ndomain,1);
ucomp = zeros(ndomain,1);
k = clmparams.k;
coef = clmparams.coef;


% compute the exact solution
for i=1:ndomain
  j=i+1;
  if j > ndomain
    j = j - ndomain;
  end
  
  uexact(i) = uexact(i) + chnk.helm2d.green(k(i),clmparams.src_in(:,j),targ(:,i));
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



% evaluate the field in the second domain at 10000 points and record time
ngr = clmparams.ngr;       % field evaluation at ngr^2 points
xylim=clmparams.xylim;  % computational domain
[xg,yg,targs,ntarg] = clm.targinit(xylim,ngr);

disp(' ')
disp(['Evaluate the field at ', num2str(ntarg), ' points'])
disp('Step 1: identify the domain for each point')
clist = clmparams.clist;
[targdomain,tid] = clm.finddomain_gui(chnkr,clmparams,targs);


tid = unique(tid);
ntid = setdiff(1:ntarg,tid);

list = cell(1,ndomain);
for i=1:ndomain
    list{i} = find(targdomain==i);
end

[sk,~,exp_mat] = clm.get_compressed_postproc_im(chnkr,clmparams);
opts = [];
tic, [eva_mats,sktarg] = clm.get_evamat_postproc_im(chnkr,clmparams,targs, ...
   targdomain,sk,eps,opts); toc






if 1==1
    disp('Step 2: evaluate the total field for point sources')

    eps0 = 1e-5;
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

alpha = pi/2;
opts_rhs = [];
opts_rhs.itype = 2;
opts_rhs.alpha = alpha;



disp(' ')
disp(['Now calculate the field when the incident wave is a plane wave'])
disp(['incident angle = ', num2str(opts_rhs.alpha)])
rhs = clm.get_rhs_gui(chnkr,clmparams,np,alpha1,alpha2,opts_rhs);
[sol] = chnk.flam.solve_2by2blk(rhs,Fskel,Fskel2,skel_struct,opts_perm);

disp('Step 4: evaluate the total field for incident plane wve')
eps0 = 1e-5;
start = tic;

%tic, [u,~] = clm.postprocess_sol_gui_fmmcorr_slower(chnkr,clmparams,targs,targdomain,eps0,sol,sk,exp_mat); toc;
tic, [u1,gradu1] = clm.postprocess_sol_gui_fmm_fds(chnkr,clmparams,targs,targdomain,eps0,sol,sk,exp_mat,eva_mats,sktarg); toc;


u = u(:);
u1 = u1(:);
start = tic; 
idomup = find(clmparams.is_inf == 1);
idomdown = find(clmparams.is_inf == -1);
for i=1:ndomain
  if ~isempty(list{i})
    if i==idomup || i==idomdown
        u(list{i}) = u(list{i}) + clm.planewavetotal_gui(k(idomup),alpha,k(idomdown),targs(:,list{i}),clmparams.is_inf(i),idomup,idomdown,coef);
        u1(list{i}) = u1(list{i}) + clm.planewavetotal_gui(k(idomup),alpha,k(idomdown),targs(:,list{i}),clmparams.is_inf(i),idomup,idomdown,coef);
    end
  end
end

dt = toc(start);


disp(['Evaluation time = ', num2str(dt), ' seconds'])
%fprintf('%5.2e s : time to evaluate the field at 10000 points\n',t1)

% plot out the field
clm.fieldplot(abs(u1),chnkr,xg,yg,xylim,ngr,fontsize)
title('Total field, incident plane wave angle = $\frac{3\pi}{4}$','Interpreter','LaTeX','FontSize',fontsize)

