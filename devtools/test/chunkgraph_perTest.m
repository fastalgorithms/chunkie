%chunkgraph_perTest0();

%function chunkgraph_perTest()
%chunkgraph_perTest: 
%
% Objectives: 
% - test geometry, region detection with plot_regions,
% chunkgraph_perinregion


%housekeeping: 
clear; close all; clc; 
vrb = true;   % set false to skip figures
if vrb
    %comp domain pts:
    Nx = 150; 
    Ny = 150; 
end

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
%{
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
    Nx = 150; Ny = 150; 
    plot_geom(cg,Nx,Ny)
end

%}

%% chunkermat RCIP: 
%
%QP problem on staircase: 
%
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg = chunkgraph_per(verts,edgesendverts,merge_idx);
dx = cg.dx; 

%src: 
src = []; src.r = [0;-0.5]; 

%geometry + source plot: 
figure; hold on; 
plot(cg); 
scatter(src.r(1),src.r(2),'ro','filled')

%computational domain: 
Nx = 100; Ny = 100; 
cell_targs = gen_comp_domain(cg,Nx,Ny); 
nper = 1; nshift = floor(nper/2); 
cg_comp = cg; 
comp_targs = cell_targs; 
dxv = [dx;0];
for s = -nshift:nshift
    if s == 0 
        continue
    end
    cg_comp = merge([cg_comp,cg + s*dxv]); 
    comp_targs.r = [comp_targs.r, cell_targs.r + s*dxv]; 
end
cell_exti = ~chunkerinterior(cg,cell_targs); 
comp_exti = repmat(cell_exti,nper,1); 

%{
figure;
plot(cg_comp); 
xc = cell_targs.r(1,:); yc = cell_targs.r(2,:); 
xcomp = comp_targs.r(1,:); ycomp = comp_targs.r(2,:); 
hold on; 
scatter(xcomp(comp_exti),ycomp(comp_exti)); 
scatter(xc(cell_exti),yc(cell_exti)); 
%}

%kappa curve: 
nkap = 60; 
dt = (2*pi/dx) / nkap; 
tkap = -pi/dx + dt*(0:nkap-1) ; 
[kap,kap_p] = kappa_curve(tkap); 
w = dt; 

%solving QP problem: 
zk = 2*dx; 
us_zk_comp = zeros(size(comp_targs.r,2),1); 
opts = []; opts.forcesmooth = true; 
parfor k = 1:nkap
    kernsp = -2*kernel('hq','sp',zk,kap(k),dx);
    rhs    = -kernsp.eval(src,cg);

    sysmat = eye(cg.npt) + chunkermat(cg,kernsp);
    sig    = sysmat\rhs;

    kerns = kernel('hq','s',zk,kap(k),dx);

    us_zk_cell = chunkerkerneval(cg,kerns,sig,cell_targs,opts); 
    us_zk_comp = us_zk_comp + w*kap_p(k) * kron(exp(1i*kap(k)*dx*(-nshift:nshift)),us_zk_cell); 
end
us = (dx/(2*pi)) * us_zk_comp; 
us(~comp_exti) = NaN; 
 

%incident wave: 
kerns = kernel('h','s',zk); 
ui = kerns.eval(src,comp_targs); 

%total field: 
u = ui + us; 

%plotting: 
plotdata = log10(abs(u)); 
plotdata(~comp_exti) = NaN; 
figure; hold on; 
xplot = reshape(xcomp,[],Nx*nper); 
yplot = reshape(ycomp,[],Ny*nper); 
plotdata = reshape(plotdata,size(xplot)); 
pcolor(xplot,yplot,plotdata)
shading interp
colorbar

%}


%test 1: 
%{
opts = [];
opts.rcip = false;
A_no_rcip = chunkermat(cg,kern,opts);

opts.rcip = true;
A_rcip = chunkermat(cg,kern,opts);

fprintf('A_no_rcip size: %d x %d\n',size(A_no_rcip,1),size(A_no_rcip,2));
fprintf('A_rcip    size: %d x %d\n',size(A_rcip,1),size(A_rcip,2));
fprintf('relative difference = %.3e\n', ...
norm(A_rcip - A_no_rcip,'fro')/norm(A_no_rcip,'fro'));
%}

%test 2: 
%{
opts = [];
opts.rcip = true;

opts.rcip_phase = false;
[A0,~,rcip0] = chunkermat(cg,kern,opts);

opts.rcip_phase = true;
[A1,~,rcip1] = chunkermat(cg,kern,opts);

for ivert = 1:numel(rcip1)

    if isempty(rcip1{ivert})
        continue
    end

    starind = rcip1{ivert}.starind;
    pedge = rcip1{ivert}.pedge;

    ngl = cg.echnks(1).k;
    ndim = 1;  % change if vector-valued kernel

    ph = repelem(pedge(:),2*ngl*ndim);

    B_expected = (ph * (1./ph).') .* A0(starind,starind);
    B_actual   = A1(starind,starind);

    relerr = norm(B_actual - B_expected,'fro')/norm(B_expected,'fro');

    fprintf('vertex %d phase block relerr = %.3e\n',ivert,relerr);
end
%}

%end

%% helpers:

function targs = gen_comp_domain(cgrph,Nx,Ny)
    verts = cgrph.verts; 
    xmin = min(verts(1,:)); xmax = max(verts(1,:)); 
    ymin = min(verts(2,:)); ymax = max(verts(2,:)); 
    x1 = linspace(xmin,xmax,Nx);
    y1 = linspace(ymin-0.5,ymax+0.5,Ny);
    [xx,yy] = meshgrid(x1,y1);
    targs = []; targs.r = [xx(:).'; yy(:).'];
end

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
    Nx = 150; Ny = 150; 
    targs = gen_comp_domain(cg,Nx,Ny); 
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

