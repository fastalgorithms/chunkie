function [isinside] = regioninside(cgrph,rgn1,rgn2)

    e2 = rgn2{1}{1}(1);
    v2 = find(cgrph.verts2edge(:,abs(e2))==1);
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
    isinside = (irgn ~=0);
end

