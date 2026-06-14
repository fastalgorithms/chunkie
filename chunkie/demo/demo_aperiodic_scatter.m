%DEMO_APERIODIC_SCATTER
%{
description: 
- set up and solve an aperiodic scattering problem on an infinite, 1D domain
- see FJ Agocs, AH Barnett paper on the method: https://arxiv.org/abs/2310.12486
- this demonstration illustrates how chunkgraph_per objects can be used to
easily set up and solve scattering problems on infinite, periodic, 1D
boundaries

author: Jonathan Shaw (jshaw6300@gmail.com)
%}

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
figure; hold on; 
plot(cg); 
scatter(src.r(1),src.r(2),'ro','filled')
hold off; 

%computational domain: 
Nxper = 7; Nyper = 1; 

%comp domain opts: 
cd_opts = []; 
cd_opts.Nx = 100; cd_opts.Ny = 150; %# pts/(period + padding)
cd_opts.pad = [0 0 0 4.5]; %[xmin xmax ymin ymax] padding outside of unit cell(s)
Nshift = floor(Nxper/2); 
[cg_comp,cell_targs,comp_targs] = gen_comp_domain(cg,Nxper,Nyper,cd_opts);
cell_eidx = ~chunkerinterior(cg,cell_targs); 
comp_eidx = repmat(cell_eidx,Nxper,1); 

%kappa curve: 
Nkap = 60; 
dt = (2*pi/dx) / Nkap; 
tkap = -pi/dx + dt*(0:Nkap-1) ; 
[kap,kap_p] = kappa_curve(tkap); 
w = dt; 

zk = 1.2; %wavenumber
us_zk_comp = zeros(size(comp_targs.r,2),1); 
opts = []; opts.forcesmooth = false; %set forcesmooth = true for speed up (bypassing near quad eval routine)

%solving integral equation + evaluating soln for each node on kappa_curve: 
parfor k = 1:Nkap
    kernsp = -2*kernel('hq','sp',zk,kap(k),dx);
    rhs    = -kernsp.eval(src,cg);

    sysmat = eye(cg.npt) + chunkermat(cg,kernsp);
    sig    = sysmat\rhs;

    kerns = kernel('hq','s',zk,kap(k),dx);

    us_zk_cell = chunkerkerneval(cg,kerns,sig,cell_targs,opts); 
    us_zk_comp = us_zk_comp + w*kap_p(k) * kron(exp(1i*kap(k)*dx*(-Nshift:Nshift)).',us_zk_cell); 
end

%scattered field:
us = (dx/(2*pi)) * us_zk_comp; 
 
%incident wave: 
kerns = kernel('h','s',zk); 
ui = kerns.eval(src,comp_targs); 

%total field: 
u = ui + us; 
u(~comp_eidx) = NaN; %set values beneath boundary to NaN

%plotting: 
Nytot = cd_opts.Ny*Nyper; Nxtot = cd_opts.Nx*Nxper; 
psize = [Nytot,Nxtot]; 
xcomp = comp_targs.r(1,:); ycomp = comp_targs.r(2,:); 
xplot = reshape(xcomp,psize); 
yplot = reshape(ycomp,psize); 

xmin = min(comp_targs.r(1,:)); xmax = max(comp_targs.r(1,:)); 
ymin = min(comp_targs.r(2,:)); ymax = max(comp_targs.r(2,:)); 
axs = [xmin xmax ymin ymax]; 

redata = reshape(real(u),psize); 
imdata = imag(u); imdata(~comp_eidx) = NaN; imdata = reshape(imdata,psize); 
abdata = reshape(abs(u),psize); 

figure; set(gcf,'theme','light'); hold on; 

subplot(1,3,1); hold on; 
pcolor(xplot,yplot,redata); shading interp; colorbar; 
scatter(src.r(1),src.r(2),'ro','filled')
plot(cg_comp,'k','linewidth',6)
axis(axs)
title('Re u','FontSize',16)
hold off; 

subplot(1,3,2); hold on; 
pcolor(xplot,yplot,imdata); shading interp; colorbar; 
plot(cg_comp,'k','linewidth',6)
scatter(src.r(1),src.r(2),'ro','filled')
axis(axs)
title('Im u','FontSize',16)
hold off; 

subplot(1,3,3); hold on; 
pcolor(xplot,yplot,abdata); shading interp; colorbar;
plot(cg_comp,'k','linewidth',6)
scatter(src.r(1),src.r(2),'ro','filled')
axis(axs)
title('|u|','FontSize',16)
hold off; 

sgtitle('Total field u')
hold off; 

%error plot for interior point source: 
%{
edata= log10(abdata);  
figure; set(gcf,'theme','light'); hold on; 
pcolor(xplot,yplot,edata); shading interp; colorbar; 
scatter(src.r(1),src.r(2),'ro','filled')
axis(axs)
plot(cg_comp,'k','linewidth',6)
title('log10(|u|)')
%}
t = toc(tstart); 

fprintf('\nElapsed time: %1.1f s\n',t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kappa curve discretization for inverse Floquet Bloch transform: 
function [kappa,kappa_p] = kappa_curve(t)
    kappa = t - 1i*sin(t); 
    kappa_p = 1 - 1i*cos(t); 
end