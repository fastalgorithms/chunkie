function poly = cell_polygon(obj,edges)
%CELL_POLYGON traces an object closed under periodic tiling by concatenating edge
%points. Used to identify periodic objects that are cut by the cell
%boundary.
%
% Syntax: poly = cell_polygon(obj,edges)
%
% Output:
%   poly - (2 x obj.npt) points on closed polygon

% author: Jonathan Shaw

    poly = []; prevend = [];
    for jj = 1:numel(edges)
        e = edges(jj);
        P = obj.echnks(abs(e)).r(1:2,:);
        if e < 0
            P = fliplr(P);
        end
        if ~isempty(prevend)
            P = P + (prevend - P(:,1));   % shift so this edge starts at prevend
        end
        poly = [poly, P];
        prevend = poly(:,end);
    end
end
