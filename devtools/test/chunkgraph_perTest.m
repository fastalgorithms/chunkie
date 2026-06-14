%chunkgraph_perTest0();

function chunkgraph_perTest()
%chunkgraph_perTest: 
%
% Objectives: 
% - test geometry, region detection with plot_regions,
% chunkgraph_perinregion
% - test chunkermat RCIP edit (shift copies locally + apply phase shift
% after)


%housekeeping: 
clear; close all; clc; 
vrb = true;   % set false to skip figures

%% geometry tests: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%object closed under x-tiling: 
%{
verts = [-0.5, -0.25, -0.5, 0.5, 0.25, 0.5; -4, -3.5, -3, -4, -3.5, -3]; 
edgesendverts = [6 5 1 2; 5 4 2 3]; 
merge_idx = {[1 4],[3 6]}; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx); 

if vrb
    plot_geom(cg,Nx,Ny)
end

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%object closed under xy-tiling: 
%{
verts = [-0.25,-0.5,-0.5,-0.25,0.25,0.5,0.5,0.25;-5,-4.5,4.5,5,5,4.5,-4.5,-5]; 
edgesendverts = [7:-2:1;8:-2:2]; 
merge_idx = {[1 4],[2 7],[3 6],[5 8]}; 
fcurve = @(t) chnk.curves.fsine(t,0.1,pi,0); 
cg = chunkgraph_per(verts,edgesendverts,merge_idx,fcurve); 

if vrb 
    plot_geom(cg,Nx,Ny)
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single staircase: 
%{
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg = chunkgraph_per(verts,edgesendverts,merge_idx);

if vrb
    plot_geom(cg,Nx,Ny)
end

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%layered staircase
%{
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg0 = chunkgraph_per(verts,edgesendverts,merge_idx);
cg  = stack_layers([cg0 + [0;-1], cg0, cg0 + [0;1]], merge_idx);

if vrb
    plot_geom(cg,Nx,Ny)
end

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%closed object with user-defined periods
%{
rng(123)
t = sort(2*pi*rand(9,1));
verts = starfish(t,5,0.3);
[~, nv] = size(verts);
edgesendverts = [1:nv; circshift(1:nv,-1)];
fchnks = []; 
cparams = []; cparams.dx = 3; cparams.dy = 3; 
merge_idx = []; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx,fchnks,cparams);

if vrb 
    plot_geom(cg,Nx,Ny); 
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%composite object: 
%
%closed object: 
verts = [0.4,-0.1,0.4;2,2.5,3]; 
[~, nv] = size(verts);
edgesendverts = [nv:-1:1; circshift(nv:-1:1,-1)];
fcurve = @(t) chnk.curves.fsine(t,0.1,2*pi,0); 
cpars = []; cpars.dx = []; cpars.dy = []; % dx, dy periodicities must be consistent with later geometries
merge_idx = []; 
cg1 = chunkgraph_per(verts,edgesendverts,merge_idx,fcurve,cpars); 

%object closed under x-tiling: 
verts = [-0.5, -0.25, -0.5, 0.5, 0.25, 0.5; -4, -3.5, -3, -4, -3.5, -3]; 
edgesendverts = [6 5 1 2; 5 4 2 3]; 
merge_idx = {[1 4],[3 6]}; 
cg2 = chunkgraph_per(verts,edgesendverts,merge_idx); 

%object closed under y-tiling: 
verts = [-0.2,0,0.2,-0.2,0,0.2;-5,-4.8,-5,5,4.8,5]; 
edgesendverts = [3 2 4 5;2 1 5 6]; 
merge_idx = {[1 4],[3 6]}; 
cg3 = chunkgraph_per(verts,edgesendverts,merge_idx); 

%object closed under xy-tiling: 
verts = [-0.25,-0.5,-0.5,-0.25,0.25,0.5,0.5,0.25;-5,-4.5,4.5,5,5,4.5,-4.5,-5]; 
edgesendverts = [7:-2:1;8:-2:2]; 
merge_idx = {[1 4],[2 7],[3 6],[5 8]}; 
fcurve = @(t) chnk.curves.fsine(t,0.1,pi,0); 
cg4 = chunkgraph_per(verts,edgesendverts,merge_idx,fcurve); 

%open periodic boundary: 
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg5 = chunkgraph_per(verts,edgesendverts,merge_idx);
cg5 = stack_layers([cg5 + [0;-1], cg5, cg5 + [0;1], cg5 + [0;4]], merge_idx);

cg = merge([cg1,cg2,cg3,cg4,cg5]); 

if vrb 
    Nx = 150; Ny = 300; 
    plot_geom(cg,Nx,Ny)
end

%}

%% chunkermat RCIP tests: 
%

%testing kappa = 0: 
%{
%chunkgraph_per: 
verts = [-0.5, 0, 0.5; -1,0,-1];
edgesendverts = [3 2; 2 1];
merge_idx = {[1 3]}; 
cgper = chunkgraph_per(verts,edgesendverts,merge_idx); 
dx = cgper.dx; 

%chunkgraph: 
verts = [-1,-0.5, 0, 0.5; 0,-1,0,-1];
edgesendverts = [4 3 2; 3 2 1];
cg = chunkgraph(verts,edgesendverts); 

zk = 1.1; 
kap = 0; 
kern = kernel('hq','s',zk,kap,dx); 

[A0,~,rcip0] = chunkermat(cg,kern); 
[A1,~,rcip1] = chunkermat(cgper,kern); 

%compare mats for cgrph vert 2, cg vert 1: 
si0 = rcip0{2}.starind; 
B0 = A0(si0,si0); 

si1 = rcip1{1}.starind; 
B1 = A1(si1,si1); 

n = size(B0,1)/2;
swap = @(B)[B(n+1:end,n+1:end) B(n+1:end,1:n); B(1:n,n+1:end) B(1:n,1:n)];
fprintf('after edge-block swap: %.3e\n', norm(B0 - swap(B1),'fro')/norm(B0,'fro')); %may have swapped indices... doesn't mean chunkermat is wrong

fprintf('no swap:   %.3e\n', norm(B0 - B1,'fro')/norm(B0,'fro'));
fprintf('swap:      %.3e\n', norm(B0 - swap(B1),'fro')/norm(B0,'fro'));
fprintf('svd match: %.3e\n', norm(sort(svd(B0))-sort(svd(B1)))/norm(svd(B0)));
%}

%full QP problem: 
%
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
src = []; src.r = [0;-0.5]; 

%geometry + source plot: 
figure; hold on; 
plot(cg); 
scatter(src.r(1),src.r(2),'ro','filled')
hold off; 

%computational domain (should choose odd Nxper): 
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

abdata = reshape(abs(u),psize); 

%error plot for interior point source: 
edata= log10(abdata);  
figure; set(gcf,'theme','light'); hold on; 
pcolor(xplot,yplot,edata); shading interp; colorbar; 
scatter(src.r(1),src.r(2),'ro','filled')
axis(axs)
plot(cg_comp,'k','linewidth',6)
title('log10(|u|)')

t = toc(tstart); 

fprintf('\nElapsed time: %1.1f s\n',t)
end

%% helpers:
function plot_geom(cg,Nx,Ny)
    %basic geometry: 
    figure; hold on; 
    
    subplot(1,3,1); hold on; 
    plot(cg); 
    axs = axis; 
    title('chunkgraph\_per geometry')
    hold off; 
    
    %plot_regions: 
    subplot(1,3,2); hold on; 
    plot_regions(cg)
    title('plot\_regions')
    
    %chunkgraph_perinregion: 
    cd_opts = []; cd_opts.Nx = Nx; cd_opts.Ny = Ny; 
    cd_opts.pad = [0 0 -1e-5 -1e-5]; 
    Nxper = 1; Nyper = 1; 
    [~,~,targs] = gen_comp_domain(cg,Nxper,Nyper,cd_opts); 
    ireg = chunkgraph_perinregion(cg,targs); 
    xx = targs.r(1,:); yy = targs.r(2,:); 
    nreg = size(cg.regions,2); 
    Legend = cell(1,nreg); 
    subplot(1,3,3); hold on; 
    for reg = 1:nreg
        reg_idx = ireg == reg; 
        scatter(xx(reg_idx),yy(reg_idx),[],ireg(reg_idx),'.'); 
        Legend{reg} = ['region',num2str(reg)]; 
    end
    leg = legend(Legend); leg.AutoUpdate = 'off'; 
    axis(axs)
    plot(cg)
    title('chunkgraph\_perinregion')
    hold off; 
    
    sgtitle('region detection')
    hold off; 
end

function [kappa,kappa_p] = kappa_curve(t)
    kappa = t - 1i*sin(t); 
    kappa_p = 1 - 1i*cos(t); 
end

