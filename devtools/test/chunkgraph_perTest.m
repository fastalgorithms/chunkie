chunkgraph_perTest0();

function chunkgraph_perTest0()
%chunkgraph_perTest: 
%
% Objectives: 
% - test geometry, region detection with plot_regions,
% chunkgraph_perinregion
% - test chunkermat RCIP edit (shift copies locally + apply phase shift
% after)

vrb = false;   % set false to skip figures

%% geometry test:
%singly-periodic, layered media: 

%composite object: 
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg = chunkgraph_per(verts,edgesendverts,merge_idx);
cg = merge([cg + [0;-1], cg, cg + [0;1], cg + [0;4]]);

if vrb 
    Nx = 150; Ny = 150; 
    plot_geom(cg,Nx,Ny)
end

%pick known pts in regions: 
Nxper = 1; Nyper = 1; 
[~,~,targs] = gen_comp_domain(cg,Nxper,Nyper); 
ireg = chunkgraph_perinregion(cg,targs); 
irknown = [0 -0.25 0 0 0; 4.5 3 0 -1 -2]; 
nr = 5; 
irtest = nan(1,nr); 
for i = 1:nr
     [~,idx] = min(vecnorm(targs.r - irknown(:,i))); 
     irtest(i) = isequal(i,ireg(idx)); 
end
assert(min(irtest)==1,'chunkgraph_perinregion test failed.')

%composite object: 
%smooth circle touching boundary:
cx = -0.45; cy = 1.5; R = 0.25;     
verts = [cx + R; cy];
edgesendverts = [1; 1];
fcurve = @(t) fcircle(t,[cx;cy],R);
cpars = []; cpars.dx = []; cpars.dy = [];
merge_idx = [];
cg0 = chunkgraph_per(verts,edgesendverts,merge_idx,fcurve,cpars);

%closed object: 
verts = [0.4,-0.1,0.4;2,2.5,3]; 
[~, nv] = size(verts);
edgesendverts = [nv:-1:1; circshift(nv:-1:1,-1)];
fcurve = @(t) chnk.curves.fsine(t,0.1,2*pi,0); 
cpars = []; cpars.dx = []; cpars.dy = [];
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
cg5 = merge([cg5 + [0;-1], cg5, cg5 + [0;1], cg5 + [0;4]]);

cg = merge([cg0,cg1,cg2,cg3,cg4,cg5]); 

if vrb 
    Nx = 150; Ny = 300; 
    plot_geom(cg,Nx,Ny)
end

%pick known pts in regions: 
Nxper = 1; Nyper = 1; 
[~,~,targs] = gen_comp_domain(cg,Nxper,Nyper); 
ireg = chunkgraph_perinregion(cg,targs); 
irknown = [0 -0.25 0 0 0 0.2 cx 0.4 0.45 0; 4.5 3 0 -1 -2 2.5 cy -3.5 -4.85 -4.9]; 
nr = 10; 
irtest = nan(1,nr); 
for i = 1:nr
     [~,idx] = min(vecnorm(targs.r - irknown(:,i))); 
     irtest(i) = isequal(i,ireg(idx)); 
end
assert(min(irtest)==1,'chunkgraph_perinregion test failed.')

%% chunkermat RCIP test: 

%Bloch phase test: 
%set up unit cell of boundary: 
verts = [-0.5, 0, 0.5; -1,0,-1];
edgesendverts = [3 2; 2 1];
merge_idx = {[1 3]}; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx); 
dx = cg.dx; 
refopts = []; refopts.nover = 1; 
cg = refine(cg,refopts); 

%src: 
src = []; src.r = [0;-0.5]; 

%computational domain (should choose odd Nxper): 
Nxper = 7; Nyper = 1; 

%comp domain opts: 
cd_opts = []; 
cd_opts.Nx = 30; cd_opts.Ny = 30; %# pts/(period + padding)
cd_opts.pad = [0 0 0 4.5]; %[xmin xmax ymin ymax] padding outside of unit cell(s)
Nshift = floor(Nxper/2); 
[cg_comp,cell_targs,comp_targs] = gen_comp_domain(cg,Nxper,Nyper,cd_opts);
cell_eidx = ~chunkerinterior(cg,cell_targs); 
comp_eidx = repmat(cell_eidx,Nxper,1); 

zk = 1.2;
opts = []; opts.forcesmooth = false; 

kap = -pi/dx; 
kernsp = -2*kernel('hq','sp',zk,kap,dx);
rhs    = -kernsp.eval(src,cg);

sysmat = eye(cg.npt) + chunkermat(cg,kernsp);
sig    = sysmat\rhs;

kerns = kernel('hq','s',zk,kap,dx);

us__cell = chunkerkerneval(cg,kerns,sig,cell_targs,opts); 
us = kron(exp(1i*kap*dx*(-Nshift:Nshift)).',us__cell); 
 
%incident wave: 
kerns = kernel('hq','s',zk,kap,dx); 
ui = kerns.eval(src,comp_targs); 

%total field: 
u = ui + us; 
meanerr = mean(abs(u(comp_eidx)),'all'); 
assert(meanerr<1e-8,'chunkermat RCIP test failed.')

if vrb

    u(~comp_eidx) = NaN; %set values beneath boundary to NaN
    Nytot = cd_opts.Ny*Nyper; Nxtot = cd_opts.Nx*Nxper; 
    psize = [Nytot,Nxtot]; 
    xcomp = comp_targs.r(1,:); ycomp = comp_targs.r(2,:); 
    xplot = reshape(xcomp,psize); 
    yplot = reshape(ycomp,psize); 
    
    xmin = min(comp_targs.r(1,:)); xmax = max(comp_targs.r(1,:)); 
    ymin = min(comp_targs.r(2,:)); ymax = max(comp_targs.r(2,:)); 
    axs = [xmin xmax ymin ymax]; 
    abdata = reshape(abs(u),psize); 
    edata= log10(abdata);  
    figure; set(gcf,'theme','light'); hold on; 
    pcolor(xplot,yplot,edata); shading interp; colorbar; 
    scatter(src.r(1),src.r(2),'ro','filled')
    axis(axs)
    plot(cg_comp,'k','linewidth',6)
    title('log10(|u|)')


end

fprintf('\nchunkgraph_per tests passed.\n')
end

%% helpers:
function [r,d,d2] = fcircle(t,c,R)
% circle of radius R centered at c = [cx;cy], parameterized over [0,1]
    ts = size(t);
    th = 2*pi*t(:);
    rx  =  c(1) + R*cos(th);          ry  =  c(2) + R*sin(th);
    dx  = -2*pi*R*sin(th);            dy  =  2*pi*R*cos(th);
    d2x = -(2*pi)^2*R*cos(th);        d2y = -(2*pi)^2*R*sin(th);
    r  = reshape([rx.';ry.'],   [2,ts]);
    d  = reshape([dx.';dy.'],   [2,ts]);
    d2 = reshape([d2x.';d2y.'], [2,ts]);
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

function [cg_comp,cell_targs,comp_targs] = gen_comp_domain(cg,Nxper,Nyper,opts)
  
    if nargin < 4
        opts = []; 
    end

    if ~isa(cg,'chunkgraph_per')
        if (~isfield(opts,'periodic')) || (~isfield(opts,'dx') && ~isfield(opts,'dy'))
        error('Object must be a chunkgraph_per or opts.periodic must be set.')
        end
    end

    dx = 0; 
    if isprop(cg,'dx') && ~isempty(cg.dx)
        dx = cg.dx;
    elseif isfield(opts,'dx')
        dx = opts.dx; 
    end

    dy = 0; 
    if isprop(cg,'dy') && ~isempty(cg.dy)
        dy = cg.dy; 
    elseif isfield(opts,'dy')
        dy = opts.dy; 
    end

    pad = zeros(1,4); 
    if isfield(opts,'pad')
        pad = opts.pad; 
    end

    if ~isfield(opts,'Nx')
        Nx = 30; 
    else
        Nx = opts.Nx; 
    end
    if ~isfield(opts,'Ny')
        Ny = 30; 
    else
        Ny = opts.Ny; 
    end
    v = cg.verts; 
    xmin = min(v(1,:)); xmax = max(v(1,:)); 
    ymin = min(v(2,:)); ymax = max(v(2,:)); 
    x1 = linspace(xmin - pad(1), xmax + pad(2), Nx); 
    y1 = linspace(ymin - pad(3), ymax + pad(4), Ny); 
    [xx_cell,yy_cell] = meshgrid(x1,y1);
    cell_targs = []; cell_targs.r = [xx_cell(:).'; yy_cell(:).']; 

    cg_comp = cg; 
    comp_targs = cell_targs; 
    Nxshift = floor(Nxper/2); 
    Nyshift = floor(Nyper/2); 
    for xshift = 1:Nxshift
        dxv = [xshift*dx;0]; 
        cg_comp = merge([cg + [-dxv],cg_comp,cg + dxv]);
        comp_targs.r = [cell_targs.r - dxv,comp_targs.r,cell_targs.r + dxv]; 
        for yshift = 1:Nyshift
            dyv = [0;yshift*dy]; 
            cg_comp = merge([cg + [-dyv],cg_comp,cg + dyv]); 
            comp_targs.r = [cell_targs.r - dyv,comp_targs.r,cell_targs.r + dyv]; 
        end
    end     
end

