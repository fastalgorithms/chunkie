%GEN_COMP_DOMAIN
%{
description: 
- generate computational domain for chunkgraph_per object

input: 
    - cg, chunkgraph or chunkgraph_per object. If object is not a
    chunkgraph_per, opts.periodic and either opts.dx or opts.dy must be set
    - Nxper, Nyper >= 1: number of object periods to plot
    - opts:
        - may contain .periodic and .dx or .dy fields for chunkgraph
        - may contain xpad or ypad to pad computational domain beyond unit
        cell
        - may contain Nx/Ny, number of points in x/y per unit cell
        
output: 
    - targs, targs.r (2 x npt) storing target points in desired
    computational domain
    - cg, chunkgraph or chunkgraph_per associated with computational domain

author: Jonathan Shaw
%}

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
        Nx = 100; 
    else
        Nx = opts.Nx; 
    end
    if ~isfield(opts,'Ny')
        Ny = 100; 
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