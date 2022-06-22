function [cgrph] = chunkgraphinit(verts,edge2verts,fchnks,prefs)
    cgrph            = chunkgraph(prefs);
    cgrph.verts      = verts;
    cgrph.edge2verts = edge2verts;
    cgrph.echnks     = chunker.empty;
    
    cparams = [];
    cparams.ta = 0;
    cparams.tb = 1;
    cparams.ifclosed = 0;
    cparams.nover = 1;
    cparams.eps = 1.0d-10;
    cparams.lvlr = 'a';
    pref = [];
    pref.nchmax = 10000;
    pref.k = 16;
    
    if (size(verts,2) ~= size(edge2verts,2))
        error('Incompatible vertex and edge sizes'); 
    end
    
    echnks = chunker.empty();
    for i=1:size(edge2verts,1)
        if (numel(fchnks)<i || isempty(fchnks{i}))
            i1 = find(edge2verts(i,:)==-1);
            i2 = find(edge2verts(i,:)==1);
            v1 = verts(:,i1);
            v2 = verts(:,i2);
            fcurve = @(t) linefunc(t,v1,v2);
            chnkr = chunkerfunc(fcurve,cparams,pref);
            echnks(i) = chnkr;
        elseif (~isempty(fchnks{i}) && isa(fchnks{i},'function_handle'))
            [vs,~,~] =fchnks{i}([0,1]);
            chnkr = chunkerfunc(fchnks{i},cparams,pref);
            vfin0 = verts(:,find(edge2verts(i,:)==-1));
            vfin1 = verts(:,find(edge2verts(i,:)== 1));
            r0 = vs(:,1);
            r1 = vfin0;
            scale = norm(vfin1-vfin0,'fro')/norm(vs(:,2)-vs(:,1),'fro');
            xdfin = vfin1(1)-vfin0(1);
            ydfin = vfin1(2)-vfin0(2);
            tfin = atan2(ydfin,xdfin);
            xdini = vs(1,2)-vs(1,1);
            ydini = vs(2,2)-vs(2,1);
            tini = atan2(ydini,xdini);
            trotat = tfin - tini;
            
            chnkr = skew(chnkr,r0,r1,trotat,scale);
            
            echnks(i) = chnkr;
        end   
    end    
    cgrph.echnks = echnks;
end

