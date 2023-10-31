function [nin] = pointinregion(cgrph,rgn,r0)
    nin = 0;
    for ii=1:numel(rgn)
        iedges = rgn{ii};
        rs = [];
        for jj=1:numel(iedges)
            rtmp = sort(cgrph.echnks(abs(iedges(jj)))).r(:,:);
            if (iedges(jj)<0)
                rtmp = fliplr(rtmp(:,:));
            end
            rs = [rs,rtmp];
        end
        [in] = inpolygon(r0(1),r0(2),rs(1,:),rs(2,:));
        if (in == 1)
            nin = nin + 1;
        end
    end
end

