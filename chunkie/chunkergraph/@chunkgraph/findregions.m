function [regions] = findregions(obj)

if (isfield(obj,'vstruc'))
    vstruc = obj.vstruc;
else
    vstruc = procverts(obj);
end    
nedge  = size(obj.edge2verts,1);
e2v    = obj.edge2verts;

edges = [1:nedge,-(1:nedge)];

regions = {};
nregions = 0;

while (numel(edges)>0)
    
    enum   = edges(1);
    edges(1) = [];
    estart = enum;
    ecycle = [enum];
    
    iv0 = find(e2v(abs(enum),:)==-sign(enum));
    ivc = find(e2v(abs(enum),:)== sign(enum));
    ifdone = false;

    while (~ifdone)
        inds = find(vstruc{ivc}{1} == abs(enum));
        irel = find(vstruc{ivc}{2}(inds) == sign(enum));
        ind  = inds(irel);
        if (ind < numel(vstruc{ivc}{1}))
            enext = vstruc{ivc}{1}(ind+1);
            esign = vstruc{ivc}{2}(ind+1);
        else
            enext = vstruc{ivc}{1}(1);
            esign = vstruc{ivc}{2}(1);
        end
        enum = -esign*enext;
        ivc = find(e2v(abs(enum),:)== sign(enum));
        if (enum == estart)
            ifdone = true;
        else
            ecycle = [ecycle,enum];
            edges(edges==enum) = [];
        end  
    end  
    
    nregions = nregions + 1;
    rcurr = {};
    rcurr{1} = ecycle;
    regions{nregions} = rcurr;
    
    
end

end

