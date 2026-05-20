function [ifin] = pointinloop(cgrph,loop,r0)
    
        iedges = loop;
        rs = [];
        for jj=1:numel(iedges)
            rtmp = sort(cgrph.echnks(abs(iedges(jj)))).r(:,:);
            if (iedges(jj)<0)
                rtmp = fliplr(rtmp(:,:));
            end
            rs = [rs,rtmp];
        end
        [ifin] = inpolygon(r0(1),r0(2),rs(1,:),rs(2,:));
end

