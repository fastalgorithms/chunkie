function [cgrph] = chunkgraphinit(verts,edge2verts,fchnks,cparams)

    warning("This method is deprecated. Use the chunkgraph constructor instead.");
    prefs = [];
    cgrph            = chunkgraph(prefs);
    cgrph.verts      = verts;
    cgrph.edge2verts = edge2verts;
    cgrph.echnks     = chunker.empty;
    
    if (nargin < 4)
        cploc = [];
        cploc.ta = 0;
        cploc.tb = 1;
        cploc.ifclosed = 0;
        cploc.nover = 1;
        cploc.eps = 1.0d-10;
        cploc.lvlr = 'a';
    else
        cploc = cparams;
        %mandatory settings
        cploc.ta = 0;
        cploc.tb = 1; 
        cploc.ifclosed = 0;
        if (~isfield(cparams,'lvlr'))
            cploc.lvlr = 'a';
        end
        if (~isfield(cparams,'eps'))
            cploc.eps = 1.0d-10;
        end
        if (~isfield(cparams,'nover'))
            cploc.nover = 1;
        end
    end
    
   
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
            fcurve = @(t) chnk.curves.linefunc(t,v1,v2);
            chnkr = chunkerfunc(fcurve,cploc,pref);
            chnkr = sort(chnkr);
            %chnkr.vert = [v1,v2];
            echnks(i) = chnkr;
        elseif (~isempty(fchnks{i}) && isa(fchnks{i},'function_handle'))
            [vs,~,~] =fchnks{i}([0,1]);
            chnkr = chunkerfunc(fchnks{i},cploc,pref);
            chnkr = sort(chnkr);
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
            chnkr = move(chnkr,r0,r1,trotat,scale);
            echnks(i) = chnkr;
        end   
    end    
    cgrph.echnks = echnks;
    cgrph.vstruc = procverts(cgrph);
    %[regions] = findregions(cgrph);
    %cgrph.regions = regions;
    
    adjmat = edge2verts'*edge2verts;
    g = graph(adjmat);
    ccomp = conncomp(g);
    
    chnkcomp = {};
    regions = {};
    
    for i=1:max(ccomp)
        inds = find(ccomp==i);
        chnkcomp{i} = inds;
        [region_comp] = findregions(cgrph,inds);
        [region_comp] = findunbounded(cgrph,region_comp);
        regions{i} = region_comp;
    end 
    
    cgrph.regions = regions;
    
end
