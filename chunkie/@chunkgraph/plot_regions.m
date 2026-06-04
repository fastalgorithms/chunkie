function plot_regions(obj,iflabel)
%PLOT_REGIONS plots regions of a chunkgraph in 2 dimensions
% All regions in the chunkgraph are plotted in a different color.
%
% Syntax: plot_regions(cgrph)
%
% Input: 
%   cgrph - chunkgraph object
%
% Optional input:
%   iflabel - integer (default 2), include a Legend showing the region
%     numbers and label edges if 2, include only the region legend if 1, 
%     no labels if 0.
%
% Output:
%   none 
%

% author: Jeremy Hoskins

if nargin < 2 || isempty(iflabel)
    iflabel = 2;
end


ifhold = ishold();

echnks =  obj.echnks;
regions = obj.regions;


nr = numel(regions);
legtext = cell(max(1,nr-1),1);

for ii=2:numel(regions)
    legtext{ii-1} = "region " + num2str(ii);
           rs = [];
        for ijk=1:numel(regions{ii}{1})
            enum = regions{ii}{1}(ijk);
            rchnk = echnks(abs(enum)).r;
            rchnk =rchnk(1:2,:);
            if (enum<0)
               rchnk = fliplr(rchnk); 
            end    
            rs = [rs, rchnk];
        end  
        if isa(obj,'chunkgraph_per') && per_region_has_unbounded(obj,regions{ii})
            % unbounded periodic region: shade the band between its two
            % bounding curves, or the half-strip beyond a single curve,
            % rather than closing an open curve into a polygon across the
            % period.
            plyrgn = per_region_poly(obj,regions{ii});
            plot(plyrgn);
            hold on
            continue
        end
        plyrgn = polyshape(rs.', 'Simplify', false);

    for jj=2:numel(regions{ii})

        for kk=2:numel(regions{ii})
        rs = [];
        for ijk=1:numel(regions{ii}{kk})
            enum = regions{ii}{kk}(ijk);
            rchnk = echnks(abs(enum)).r;
            rchnk =rchnk(1:2,:);
            if (enum<0)
               rchnk = fliplr(rchnk); 
            end    
            rs = [rs,rchnk];
        end
        
        plyrgnsub = polyshape(rs', 'Simplify', false);
        plyrgn = subtract(plyrgn, plyrgnsub);
        end  
        
    end    
    plot(plyrgn);
    hold on
end    

if nr > 1 && iflabel > 0
    legend(legtext,'AutoUpdate','off');
end

plot(obj,'k-');

if iflabel > 1
    rmin = min(obj.r(:,:), [], 2);
    diam = max(obj.r(:,:)-rmin, [], 2);
    rlbl = rmin-0.1*diam;
    % for an unbounded periodic geometry region 1 is the upper half-space,
    % so label it at the top rather than the bottom-left corner.
    if isa(obj,'chunkgraph_per') && per_any_unbounded(obj)
        rlbl(2) = max(obj.r(2,:)) + 0.2*diam(2);
    end
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

function tf = per_edges_unbounded(obj,edgelist)
%PER_EDGES_UNBOUNDED true if the boundary edge list has nonzero net
% displacement (a curve that only closes through periodic identification).
    tf = norm(loop_displacement(obj,edgelist)) > 1e-10;
end


function tf = per_any_unbounded(obj)
    tf = false;
    for ii = 1:numel(obj.regions)
        comp = obj.regions{ii};
        if iscell(comp) && ~isempty(comp) && ~isempty(comp{1})
            if per_edges_unbounded(obj,comp{1})
                tf = true; return
            end
        end
    end
end


function ply = per_halfstrip(obj,edgelist,rs)
%PER_HALFSTRIP polygon covering the half-strip on the normal side of a
% single unit-cell curve, bounded by a horizontal line well outside the
% curve so it fills the visible axes after clipping. Used for the top and
% bottom (unbounded) regions.
    ny = loop_normal_y(obj,edgelist);

    xs = rs(1,:); ys = rs(2,:);
    pad = 5*max(max(ys)-min(ys), max(xs)-min(xs));
    if ny >= 0
        ylev = max(ys) + pad;   % normal points up -> fill upward
    else
        ylev = min(ys) - pad;   % normal points down -> fill downward
    end
    pts = [rs, [xs(end); ylev], [xs(1); ylev]];
    ply = polyshape(pts.', 'Simplify', false);
end


function tf = per_region_has_unbounded(obj,comps)
%PER_REGION_HAS_UNBOUNDED true if any component of a region is an unbounded
% periodic curve.
    tf = false;
    if ~iscell(comps)
        return
    end
    for c = 1:numel(comps)
        if ~isempty(comps{c}) && per_edges_unbounded(obj,comps{c})
            tf = true; return
        end
    end
end


function ply = per_region_poly(obj,comps)
%PER_REGION_POLY polygon for an unbounded periodic region: a band between
% the upper and lower bounding curves when there are two, or a half-strip
% beyond a single curve (top/bottom regions). Closed components, if any,
% are ignored for now (relevant only to composite geometries).
    curves = {}; meany = [];
    for c = 1:numel(comps)
        if per_edges_unbounded(obj,comps{c})
            curves{end+1} = comps{c};
            meany(end+1)  = curve_mean_y(obj,comps{c});
        end
    end
    if numel(curves) >= 2
        [~,iu] = max(meany);
        [~,il] = min(meany);
        ply = per_band(obj,curves{iu},curves{il});
    else
        rs = curve_points(obj,curves{1});
        ply = per_halfstrip(obj,curves{1},rs);
    end
end


function ply = per_band(obj,e_up,e_lo)
%PER_BAND polygon for the strip between an upper and a lower curve, each a
% graph over the periodic axis (sorted by x so the band is a simple
% polygon over one period).
    ru = curve_points(obj,e_up);
    rl = curve_points(obj,e_lo);
    [~,iu] = sort(ru(1,:)); ru = ru(:,iu);
    [~,il] = sort(rl(1,:)); rl = rl(:,il);
    pts = [ru, fliplr(rl)];
    ply = polyshape(pts.', 'Simplify', false);
end


function rs = curve_points(obj,edges)
%CURVE_POINTS concatenated boundary points of an edge list.
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


function y = curve_mean_y(obj,edges)
%CURVE_MEAN_Y mean y-coordinate of the points making up an edge list.
    s = 0; cnt = 0;
    for jj = 1:numel(edges)
        rr = obj.echnks(abs(edges(jj))).r(:,:);
        s = s + sum(rr(2,:));
        cnt = cnt + size(rr,2);
    end
    y = s/cnt;
end
