%DEMO_APERIODIC_SCATTER
% 
% description: 
% - set up and solve an aperiodic scattering problem on an infinite, 1D domain
% - see FJ Agocs, AH Barnett paper on the method: https://arxiv.org/abs/2310.12486
% - this demonstration illustrates how chunkgraph_per objects can be used to
% easily set up and solve scattering problems on infinite, periodic, 1D
% boundaries

%set up unit cell of boundary: 
tstart = tic; 
verts = [-0.5, 0, 0.5; -1,0,-1];
edgesendverts = [3 2; 2 1];
merge_idx = {[1 3]}; %merge_idx = cell array, stores vertices to be periodically identified
cg = chunkgraph_per(verts,edgesendverts,merge_idx); 
dx = cg.dx; 
refopts = []; refopts.nover = 1; 
cg = refine(cg,refopts); 

%src: 
src = []; src.r = [-.25;1]; 

%geometry + source plot: 
figure(1); hold on; 
plot(cg); 
scatter(src.r(1),src.r(2),'ro','filled')
hold off; 

%computational domain: 
Nx = 100; Ny = 100; %pts per x,y period
v = cg.verts; 
xmin = min(v(1,:)); xmax = max(v(1,:)); 
ymin = min(v(2,:)); ymax = max(v(2,:)); 
pad = [0 0 0 4.5]; %[xmin xmax ymin ymax] padding outside of unit cell(s)
x1 = linspace(xmin - pad(1), xmax + pad(2), Nx); 
y1 = linspace(ymin - pad(3), ymax + pad(4), Ny); 
[xx_cell,yy_cell] = meshgrid(x1,y1);
cell_targs = []; cell_targs.r = [xx_cell(:).'; yy_cell(:).']; 
cell_in = ~chunkerinterior(cg,cell_targs); 
cell_eval_targs.r = cell_targs.r(:,cell_in); 

cg_comp = cg; 
comp_targs = cell_targs; 
for xshift = 1:5
    dxv = [xshift*dx;0]; 
    cg_comp = merge([cg + [-dxv],cg_comp,cg + dxv]);
    comp_targs.r = [cell_targs.r - dxv,comp_targs.r,cell_targs.r + dxv]; 
end     
comp_in = repmat(cell_in,1,11); 
comp_eval_targs = []; comp_eval_targs.r = comp_targs.r(:,comp_in); 

%kappa curve: 
Nkap = 60; 
dt = (2*pi/dx) / Nkap; 
tkap = -pi/dx + dt*(0:Nkap-1) ; 
[kap,kap_p] = kappa_curve(tkap); 
w = dt; 

zk = 1.2; %wavenumber
us_zk_comp = zeros(size(comp_eval_targs.r,2),1); 
opts = []; opts.forcesmooth = false; %set forcesmooth = true for speed up (bypassing near quad eval routine)

%solving integral equation + evaluating soln for each node on kappa_curve: 
parfor k = 1:Nkap
    kernsp = -2*kernel('hq','sp',zk,kap(k),dx);
    rhs    = -kernsp.eval(src,cg);

    sysmat = eye(cg.npt) + chunkermat(cg,kernsp);
    sig    = sysmat\rhs;

    kerns = kernel('hq','s',zk,kap(k),dx);

    us_zk_cell = chunkerkerneval(cg,kerns,sig,cell_eval_targs,opts); 
    us_zk_comp = us_zk_comp + w*kap_p(k) * kron(exp(1i*kap(k)*dx*(-5:5)).',us_zk_cell); 
end

%scattered field:
us = nan(numel(comp_targs.r(1,:)),1); 
us(comp_in) = (dx/(2*pi)) * us_zk_comp; 
 
%incident wave: 
kerns = kernel('h','s',zk); 
ui = kerns.eval(src,comp_targs); 

%total field: 
u = ui + us; 

%plotting: 
psize = [Ny,11*Nx]; 
xcomp = comp_targs.r(1,:); ycomp = comp_targs.r(2,:); 
xplot = reshape(xcomp,psize); 
yplot = reshape(ycomp,psize); 

xmin = min(comp_targs.r(1,:)); xmax = max(comp_targs.r(1,:)); 
ymin = min(comp_targs.r(2,:)); ymax = max(comp_targs.r(2,:)); 
axs = [xmin xmax ymin ymax]; 

redata = reshape(real(u),psize); 
imdata = imag(u); imdata(~comp_in) = NaN; imdata = reshape(imdata,psize); 
abdata = reshape(abs(u),psize); 

figure(2); set(gcf,'theme','light'); hold on; 
pcolor(xplot,yplot,imdata); shading interp; colorbar; 
plot(cg_comp,'k','linewidth',6)
scatter(src.r(1),src.r(2),'ro','filled')
axis(axs)
axis equal 
title('Im u','FontSize',16)

t = toc(tstart); 

fprintf('\nElapsed time: %1.1f s\n',t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kappa curve discretization for inverse Floquet Bloch transform: 
function [kappa,kappa_p] = kappa_curve(t)
    kappa = t - 1i*sin(t); 
    kappa_p = 1 - 1i*cos(t); 
end