function m = loop_max_jump(obj,edges)
%LOOP_MAX_JUMP largest gap between the end of one edge and the start of the
% next around a loop (with wraparound). ~0 for a geometrically continuous
% loop; ~a period for a loop that closes through periodic identification
% (period jumps). 
%
% Syntax: m = loop_max_jump(obj,edges)

    ne = numel(edges);
    if ne == 0
        m = 0; return
    end
    starts = zeros(2,ne); ends = zeros(2,ne);
    for jj = 1:ne
        e = edges(jj);
        ech = obj.echnks(abs(e));
        [r1,~] = chunkends(ech,1);
        [r2,~] = chunkends(ech,ech.nch);
        if e > 0
            starts(:,jj) = r1(:,1); ends(:,jj) = r2(:,2);
        else
            starts(:,jj) = r2(:,2); ends(:,jj) = r1(:,1);
        end
    end
    m = 0;
    for jj = 1:ne
        nx = mod(jj,ne) + 1;
        g = norm(starts(:,nx) - ends(:,jj));
        if g > m
            m = g;
        end
    end
end
