function ny = loop_normal_y(obj,edges)
%LOOP_NORMAL_Y average y-component of the orientation-adjusted normal along
% an edge list. Used to orient a curve so its normal points "up" (toward
% the upper half-space). A negative-signed edge flips the stored normal.
%
% Syntax: ny = loop_normal_y(obj,edges)
%
% Input:
%   obj   - a chunkgraph_per object
%   edges - signed edge-index list
%
% Output:
%   ny - mean y-component of the (orientation-adjusted) normal
%
% author: Jonathan Shaw

    s = 0; cnt = 0;
    for jj = 1:numel(edges)
        e = edges(jj);
        nn = obj.echnks(abs(e)).n(:,:);
        nyj = nn(2,:);
        if e < 0
            nyj = -nyj;
        end
        s = s + sum(nyj);
        cnt = cnt + numel(nyj);
    end
    ny = s/cnt;
end
