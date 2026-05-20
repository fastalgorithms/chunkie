function [isinside] = loopinside(cgrph,lpi,lpo)
    
    ei = lpi;
    eo = lpo;
    
    if (numel(intersect(ei,eo)) ~= 0)
        isinside = false;
        return
    end

    vi = cgrph.edgesendverts(2,abs(ei));
    vo = cgrph.edgesendverts(2,abs(eo));

    if (numel(intersect(vi,vo)) ~= 0)
        isinside = false;
        return
    end

    eo = lpo(1);
    vo = cgrph.edgesendverts(2,abs(eo));
    vo = cgrph.verts(:,vo);
    irgn = 0;
    isinside = pointinloop(cgrph,lpi,vo);

end

