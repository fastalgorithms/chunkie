function [unbnd, bnd] = classify_loops(obj, loops)
%CLASSIFY_LOOPS split candidate boundary loops into unbounded periodic
% curves and bounded/closed objects.
%
% Syntax:
%   [unbnd, bnd] = classify_loops(obj, loops)
%
% Input:
%   obj   - a chunkgraph_per object
%   loops - cell array of candidate signed edge lists. Duplicates and both
%           orientations of a loop are allowed; they are reduced internally.
%
% Output:
%   unbnd - cell array of unbounded periodic curves, each oriented so its
%           normal points up (+y), sorted by descending mean y
%   bnd   - cell array of closed objects, each oriented counter-clockwise
%           (positive signed area), sorted by descending mean y
%
% author: Jonathan Shaw

    dtol = 1e-10;

    loops = unique_loops(loops);

    unbnd = {}; cmy = [];
    bnd   = {}; smy = [];

    for k = 1:numel(loops)
        e = loops{k};
        if isempty(e)
            continue
        end

        if norm(loop_displacement(obj,e)) > dtol
            if loop_normal_y(obj,e) < 0
                e = -fliplr(e);
            end
            unbnd{end+1} = e;                  
            cmy(end+1)   = loop_mean_y(obj,e); 
        else
            poly = cell_polygon(obj,e);
            x = poly(1,:); y = poly(2,:);
            A = 0.5*sum(x.*y([2:end 1]) - x([2:end 1]).*y); %signed area
            if A < 0
                e = -fliplr(e);
            end
            bnd{end+1} = e;                    
            smy(end+1) = mean(y); 
        end
    end

    [~,oc] = sort(cmy,'descend');
    unbnd = unbnd(oc);

    [~,os] = sort(smy,'descend');
    bnd = bnd(os);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loops = unique_loops(loops_in)
%UNIQUE_LOOPS drop empty loops and the duplicate reverse orientation loops 
    loops = {};
    keys  = {};
    for k = 1:numel(loops_in)
        e = loops_in{k};
        if isempty(e)
            continue
        end
        key = sort(abs(e));
        isnew = true;
        for j = 1:numel(keys)
            if isequal(keys{j},key)
                isnew = false;
                break
            end
        end
        if isnew
            keys{end+1}  = key;
            loops{end+1} = e;   
        end
    end
end

