function flagslf = flagself(srcs, targs, tol)
    % identify sources and targets pairs that are within tol (1e-14);

    % todo: make this more robust by searching for all pairs that are close
    % not just the closest
    if nargin < 3
        tol = 1e-14;
    end

    flagslf = [];

    randangle = 2*pi*rand();
    randrot = [cos(randangle), -sin(randangle); ...
        sin(randangle), cos(randangle)];
    
    srcrot  = randrot*srcs(:,:);
    targrot = randrot*targs(:,:);

    [srcsortx, jds] = sort(srcrot(1,:));
    [targsortx, ids] = sort(targrot(1,:));

    binids = cell(length(srcsortx),1);

     idcheck = [1;2];
     for j = 1:length(srcsortx)
        while (idcheck(2) < length(targsortx) && ...
            targsortx(idcheck(2)) < srcsortx(j))
            idcheck = idcheck+1;
        end
        [d, kd] = min(abs(targsortx(idcheck) - srcsortx(j)));
        idclose = idcheck(kd);
        if d < tol
            binids{jds(j)} = [binids{jds(j)}, ids(idclose)];
        end
     end  

    [srcsorty, jds] = sort(srcrot(2,:));
    [targsorty, ids] = sort(targrot(2,:));
    
    idcheck = [1;2];
    for j = 1:length(srcsorty)
        while (idcheck(2) < length(targsorty) && ...
                targsorty(idcheck(2)) < srcsorty(j))
            idcheck = idcheck+1;
        end
        [d, kd] = min(abs(targsorty(idcheck) - srcsorty(j)));
        idclose = idcheck(kd);
        if d < tol
            if ismember(ids(idclose), binids{jds(j)})
                flagslf = [flagslf, [jds(j);ids(idclose)]];
            end
        end
     end  
     if ~isempty(flagslf)
         [~, ids] = sort(flagslf(1,:));
         flagslf = flagslf(:,ids);
     end
end