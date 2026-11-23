%DEMO_ELAST2D_SQUARE_DIRICHLET
%
% Interior displacement (Dirichlet) problem of 2D linear elasticity on
% the square [-1,1]^2, with displacement (1,1) prescribed on the top edge
% and zero displacement on the other three edges.
%
% Representation: u = Dalt[sigma], the alternative (principal value free)
% elastic double layer of kernel.elast2d('dalt',...). Its interior limit
% is c*sigma + Dalt sigma with c = (lam+3*mu)/(2*(lam+2*mu)) (as in
% elastickernelsTest), and on a straight edge the kernel vanishes
% identically. RCIP in chunkermat assumes a system of the form I + K, so
% we solve
%
%       (I + Dalt/c) sigma = g/c.
%
% Note: the traction with respect to the source of the string field of
% kernel.elast2d_string is a valid (Navier) double layer, but its jump is
% C = -((lam+2mu) tau tau^T + mu n n^T)/(lam+mu) and on a straight edge its
% principal value part is M/(t-s) with det(C -+ i*pi*M) = 0 identically,
% so the resulting integral equation is not Fredholm. Dalt does not have
% this problem.
%
% The script first checks the solver against an exact solution (a sum of
% Kelvinlets placed outside the square), then solves the problem above
% and plots u_1, u_2 and |u|.

% BEGIN DEMO ELAST2D SQUARE DIRICHLET
%% Lame parameters

lam = 1.5;
mu = 2.1;
c = (lam+3*mu)/(2*(lam+2*mu));

%% geometry: square as a chunkgraph (counterclockwise, outward normals)

verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];   % edge 3 (vertex 3 -> 4) is the top
cparams = []; cparams.nchmin = 8;
pref = []; pref.k = 16;
cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);
npt = cgrph.npt;
rb = cgrph.r(:,:);

%% kernels and system matrix (with RCIP at the corners)

kernd = kernel.elast2d('dalt',lam,mu);

opts = []; opts.rcip = true;
start = tic; sysmat = chunkermat(cgrph,(1/c)*kernd,opts); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1);
sys = eye(2*npt) + sysmat;

%% check against an exact solution: Kelvinlets outside the square

rng(8675309,'twister');
ns = 6;
ts = 2*pi*rand(1,ns);
src = []; src.r = 2.5*[cos(ts); sin(ts)];
f = randn(2*ns,1);
kerns = kernel.elast2d('s',lam,mu);

bdry = []; bdry.r = rb;
gex = kerns.eval(src,bdry)*f;
sigex = gmres(sys,gex/c,[],1e-13,200);

ntest = 200;
tt = 2*pi*rand(1,ntest); rr = 0.95*sqrt(rand(1,ntest));
targtest = [rr.*cos(tt); rr.*sin(tt)];
tinfo = []; tinfo.r = targtest;
uex = kerns.eval(src,tinfo)*f;
utest = chunkerkerneval(cgrph,kernd,sigex,targtest);
fprintf('Kelvinlet check: relative error %5.2e\n', ...
    norm(utest(:)-uex)/norm(uex));

%% boundary data: (1,1) on the top edge, zero elsewhere

itop = abs(rb(2,:) - 1) < 1e-12;
g = zeros(2,npt);
g(:,itop) = 1;
g = g(:);

start = tic; [sigma,~,relres,iter] = gmres(sys,g/c,[],1e-13,200);
t1 = toc(start);
fprintf('%5.2e s : time for gmres, %d iters, relres %5.2e\n',t1, ...
    iter(end),relres);

%% evaluate on a grid inside the square

nplot = 200;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,targs);
targin = targs(:,in);
nt = size(targin,2);

start = tic; u = chunkerkerneval(cgrph,kernd,sigma,targin);
t1 = toc(start);
fprintf('%5.2e s : time to evaluate at %d targets\n',t1,nt);
u = reshape(u,2,nt);

%% plot

fplot = {u(1,:), u(2,:), vecnorm(u)};
tlabs = {'$u_1$', '$u_2$', '$|u|$'};

figure(1); clf
tiledlayout(1,3,'TileSpacing','compact');
for ip = 1:3
    nexttile
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz);
    set(h,'EdgeColor','none')
    hold on
    plot(cgrph,'k','LineWidth',2)
    plot([-1 1],[1 1],'r','LineWidth',3)
    if ip == 3
        u1 = nan(size(xx)); u1(in) = u(1,:);
        u2 = nan(size(xx)); u2(in) = u(2,:);
        isk = 1:12:nplot;
        quiver(xx(isk,isk),yy(isk,isk),u1(isk,isk),u2(isk,isk),'w');
    end
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
sgtitle('displacement $(1,1)$ on top (red), zero elsewhere', ...
    'Interpreter','latex')
% END DEMO ELAST2D SQUARE DIRICHLET
