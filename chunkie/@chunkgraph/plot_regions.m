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
        if isa(obj,'chunkgraph_per') && per_edges_unbounded(obj,regions{ii}{1})
            % unbounded periodic region: shade the half-strip on the normal
            % side of the unit-cell curve instead of closing the open curve
            % into a (self-intersecting) polygon across the period.
            plyrgn = per_halfstrip(obj,regions{ii}{1},rs);
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
    tf = norm(per_net_disp(obj,edgelist)) > 1e-10;
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


function d = per_net_disp(obj,edgelist)
    d = [0;0];
    for ie = 1:numel(edgelist)
        e = edgelist(ie);
        ech = obj.echnks(abs(e));
        [r1,~] = chunkends(ech,1);
        [r2,~] = chunkends(ech,ech.nch);
        if e > 0
            d = d + (r2(:,2) - r1(:,1));
        else
            d = d + (r1(:,1) - r2(:,2));
        end
    end
end


function ply = per_halfstrip(obj,edgelist,rs)
%PER_HALFSTRIP polygon covering the half-strip on the normal side of the
% unit-cell curve, bounded above/below by a horizontal line well outside
% the curve so it fills the visible axes after clipping.
    s = 0; cnt = 0;
    for ie = 1:numel(edgelist)
        e = edgelist(ie);
        nn = obj.echnks(abs(e)).n(:,:);
        nyj = nn(2,:);
        if e < 0
            nyj = -nyj;
        end
        s = s + sum(nyj);
        cnt = cnt + numel(nyj);
    end
    ny = s/cnt;

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
