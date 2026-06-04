function poly = cell_polygon(obj,edges)
%CELL_POLYGON contiguous closed polygon traced by a periodic cell that is
% closed under tiling (e.g. the diamond of case (b)). The period jumps
% between consecutive edges are removed by shifting each edge so it starts
% where the previous one ended, yielding a single closed polygon. Only
% meaningful for loops with ~zero net displacement.
%
% Syntax: poly = cell_polygon(obj,edges)
%
% Output:
%   poly - 2 x M array of polygon vertices

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
