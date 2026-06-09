function plot_regions(obj,iflabel)
%PLOT_REGIONS plots regions of a chunkgraph_per in 2 dimensions. 
%
% Syntax: plot_regions_per(cgrph_per)
%
% Input:
%   cgrph_per - chunkgraph_per object
%
% Optional input:
%   iflabel - integer (default 2), include a legend showing region numbers
%     and label edges if 2, include only the region legend if 1, no labels
%     if 0.
%
% Output:
%   none
% authors: Jeremy Hoskins, Jonathan Shaw

if nargin < 2 || isempty(iflabel)
    iflabel = 2;
end

msg = "plot_regions_per: input 1 must be chunkgraph_per";
assert(class(obj) == "chunkgraph_per",msg);

ifhold = ishold();

regions = obj.regions;
nr = numel(regions);
legtext = cell(max(1,nr-1),1);

for ii = 2:nr
    legtext{ii-1} = "region " + num2str(ii);

    comps = regions{ii};
    if isempty(comps)
        continue
    end

    if per_reg_has_unb(obj,comps)
        plyrgn = per_reg_poly(obj,comps);
    elseif cl_und_tiling(obj,comps)
        plyrgn = per_cl_poly(obj,comps);
    else
        plyrgn = cl_poly(obj,comps);
    end

    plot(plyrgn);
    hold on
end

if nr > 1 && iflabel > 0
    legend(legtext,'AutoUpdate','off');
end

plot(obj,'k-');

if iflabel > 1
    xmin = min(obj.r(1,:)); xmax = max(obj.r(1,:)); 
    bf = 0.1; 
    xbuff = bf*(xmax-xmin); 
    ymin = min(obj.r(2,:)); ymax = max(obj.r(2,:)); 
    ybuff = bf*(ymax-ymin); 
    rlbl = [xmin-xbuff,ymax+ybuff]; 
    text(rlbl(1), rlbl(2), 'region 1');

    nedge = length(obj.echnks);
    for j = 1:nedge
        wts = obj.echnks(j).wts;
        l1 = sum(wts(:));
        l2 = cumsum(wts(:));
        ind = find(l2 > l1/2,true,'first');
        r = obj.echnks(j).r(:,ind);
        text(r(1),r(2),"edge "+num2str(j),'HorizontalAlignment','center')
    end
end

hold off

if ifhold
    hold on
end

end

function ply = cl_poly(obj,comps)
%CL_PLY polyshape for a bounded ordinary region. The first
% component is the outer boundary and later components are holes.
    ply = polyshape();
    if isempty(comps) || isempty(comps{1})
        return
    end

    rs = curve_points(obj,comps{1});
    ply = polyshape(rs.', 'Simplify', false);

    for kk = 2:numel(comps)
        if isempty(comps{kk})
            continue
        end
        rs = curve_points(obj,comps{kk});
        plysub = polyshape(rs.', 'Simplify', false);
        ply = subtract(ply,plysub);
    end
end

function eu = per_edge_unb(obj,edgelist)
%PER_EDGE_UNB true if the boundary edge list has nonzero net
% displacement, i.e. it closes only through a periodic translation.
    dtol = 1e-10; 
    eu = norm(loop_displacement(obj,edgelist)) > dtol;
end

function hb = per_reg_has_unb(obj,comps)
%PER_REG_HAS_UNB true if any component of a region is an unbounded
% periodic curve.
    hb = false;
    if ~iscell(comps)
        return
    end
    for c = 1:numel(comps)
        if ~isempty(comps{c}) && per_edge_unb(obj,comps{c})
            hb = true;
            return
        end
    end
end

function ply = per_reg_poly(obj,comps)
%PER_REG_POLY polygon for an unbounded periodic region: a band between
% two bounding curves, or a half-strip beyond a single curve. Closed
% components in the same region are ignored here because their interiors are
% plotted as their own overriding regions.
    curves = {};
    meany = [];
    for c = 1:numel(comps)
        if ~isempty(comps{c}) && per_edge_unb(obj,comps{c})
            curves{end+1} = comps{c}; 
            meany(end+1)  = loop_mean_y(obj,comps{c}); 
        end
    end

    if isempty(curves)
        ply = polyshape();
    elseif numel(curves) >= 2
        [~,iu] = max(meany);
        [~,il] = min(meany);
        ply = per_band(obj,curves{iu},curves{il});
    else
        rs = curve_points(obj,curves{1});
        ply = per_halfstrip(obj,curves{1},rs);
    end
end

function ply = per_halfstrip(obj,edgelist,rs)
%PER_HALFSTRIP polygon covering the half-strip on the normal side of a
% single periodic curve. Used for the top and bottom layered regions.
    ny = loop_normal_y(obj,edgelist);
    xs = rs(1,:);
    ys = rs(2,:);

    ptol = 1e-13; 
    if obj.dy>ptol
        if ny >= 0
            pad = max(obj.r(2,:)) - max(ys); 
        else
            pad = min(ys) - min(obj.r(2,:)); 
        end
    else
        pad = 5*max(max(ys)-min(ys), max(xs)-min(xs));
        if pad == 0
            pad = 1;
        end
    end

    if ny >= 0
        ylev = max(ys) + pad;   
    else
        ylev = min(ys) - pad;  
    end

    pts = [rs, [xs(end); ylev], [xs(1); ylev]];
    ply = polyshape(pts.', 'Simplify', false);
end

function ply = per_band(obj,e_up,e_lo)
%PER_BAND polygon for the strip between an upper and a lower curve.
    ru = curve_points(obj,e_up);
    rl = curve_points(obj,e_lo);
    [~,iu] = sort(ru(1,:));
    [~,il] = sort(rl(1,:));
    ru = ru(:,iu); rl = rl(:,il); 
    pts = [ru, fliplr(rl)];
    ply = polyshape(pts.', 'Simplify', false);
end

function clut = cl_und_tiling(obj,comps)
%CL_UND_TILING true if any component closes only through
% periodic identification: zero net displacement but nontrivial period jump.
    dtol = 1e-10; 
    jtol = 1e-6; 
    clut = false;
    if ~iscell(comps)
        return
    end
    for c = 1:numel(comps)
        e = comps{c};
        if ~isempty(e) && norm(loop_displacement(obj,e)) < dtol ...
                && loop_max_jump(obj,e) > jtol
            clut = true;
            return
        end
    end
end

function ply = per_cl_poly(obj,comps)
%PER_CL_POLY closed-under-tiling cell(s), shown over one unit cell.
% Complete cells can straddle a period boundary, so period translates are
% added and clipped to the displayed unit cell.
    
    ptol = 1e-13; 
    dtol = 1e-10; jtol = 1e-6; 
    ply = polyshape(); 
    if obj.dx>ptol && obj.dy>ptol
        for c = 1:numel(comps)
            e = comps{c}; 
            if norm(loop_displacement(obj,e)) < dtol && loop_max_jump(obj,e) > jtol
                base = cell_polygon(obj,e); 
                for j = -1:1
                    for k = -1:1 %construct 2 surrounding copies for closure:
                        ply = union(ply, polyshape((base + [j*obj.dx;k*obj.dy]).', 'Simplify', false));
                    end
                end
            end
        end
    elseif obj.dx>ptol
        for c = 1:numel(comps)
            e = comps{c}; 
            if norm(loop_displacement(obj,e)) < dtol && loop_max_jump(obj,e) > jtol
                base = cell_polygon(obj,e); 
                for j = -1:1 %construct 2 surrounding copies for closure:
                    ply = union(ply, polyshape((base + [j*obj.dx;0]).', 'Simplify', false));
                end
            end
        end
    else
        for c = 1:numel(comps)
            e = comps{c}; 
            if norm(loop_displacement(obj,e)) < dtol && loop_max_jump(obj,e) > jtol
                base = cell_polygon(obj,e); 
                for k = -1:1 %construct 2 surrounding copies for closure:
                    ply = union(ply, polyshape((base + [0;k*obj.dy]).', 'Simplify', false));
                end
            end
        end
    end

    %clipping to unit cell: 
    xmin = min(obj.r(1,:)); xmax = max(obj.r(1,:)); 
    ymin = min(obj.r(2,:)); ymax = max(obj.r(2,:)); 
    uc = polyshape([xmin xmax xmax xmin],[ymin ymin ymax ymax]); 
    ply = intersect(ply,uc);
end

function rs = curve_points(obj,edges)
%CURVE_POINTS concatenated boundary points of a signed edge list.
    rs = [];
    for ijk = 1:numel(edges)
        enum = edges(ijk);
        rchnk = obj.echnks(abs(enum)).r;
        rchnk = rchnk(1:2,:);
        if enum < 0
            rchnk = fliplr(rchnk);
        end
        rs = [rs, rchnk]; 
    end
end

