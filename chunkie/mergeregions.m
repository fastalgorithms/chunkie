function [rgnout] = mergeregions(cgrph,rgn1,rgn2)

    rgnout = {};
    
    e2 = rgn2{1}{1}(1);
    v2 = find(cgrph.edge2verts(:,abs(e2))==1);
    v2 = cgrph.verts(:,v2);

    irgn = 0;
    for ii=2:numel(rgn1)
        nin = pointinregion(cgrph,rgn1{ii},v2);
        disp("first inclusion:")
        nin
        if (nin > 0 && mod(nin,2)==1)
            if (irgn ~= 0)
                disp("Warning: an unsupported geometry error has occurred");
            end
            irgn = ii;
        end
    end

    if (irgn ~= 0)
       rgnout = rgn1;
       rgnout = [rgnout,rgn2(2:end)];
       rgnout{irgn} = [rgnout{irgn},rgn2{1}];
    end
    
    irgn1 = irgn;
    
    if (irgn == 0)
        e1 = rgn1{1}{1}(1);
        v1 = find(cgrph.edge2verts(:,abs(e1))==1);
        v1 = cgrph.verts(:,v1);

        irgn = 0;
        for ii=2:numel(rgn2)
            nin = pointinregion(cgrph,rgn2{ii},v1);
            disp("second inclusion:")
            nin
            if (nin > 0 && mod(nin,2)==1)
                if (irgn ~= 0)
                    disp("Warning: an unsupported geometry error has occurred");
                end
                irgn = ii;
            end
        end
        
        if (irgn ~= 0)
            disp("here")
            irgn
            rgnout = rgn2;
            rgnout = [rgnout,rgn1(2:end)];
            rgnout{irgn} = [rgnout{irgn},rgn1{1}];
        end
    end

    if (irgn1==0 && irgn==0)
        disp("here also")
        rgnout = rgn1;
        rgnout = [rgnout,rgn2(2:end)];
        rgnout{1} = [rgnout{1},rgn2{1}];
    end    
    
end

