function d = loop_displacement(obj,edges)
%LOOP_DISPLACEMENT net displacement around a loop, from edge chunker
% endpoints. ~0 for a genuinely closed loop; a nonzero lattice vector
% (+/-dx,0)/(0,+/-dy) for a curve that only closes through the periodic
% identification (e.g. a staircase unit cell).
%
% Syntax: d = loop_displacement(obj,edges)
%
% Input:
%   obj   - a chunkgraph object
%   edges - signed edge-index list describing a loop
%
% Output:
%   d - 2x1 net displacement vector

    d = [0;0];
    for jj = 1:numel(edges)
        e = edges(jj);
        ech = obj.echnks(abs(e));
        [r1,~] = chunkends(ech,1);
        [r2,~] = chunkends(ech,ech.nch);
        rstart = r1(:,1);
        rend   = r2(:,2);
        if e > 0
            d = d + (rend - rstart);
        else
            d = d + (rstart - rend);
        end
    end
end
