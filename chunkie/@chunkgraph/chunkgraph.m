classdef chunkgraph
%CHUNKGRAPH chunk graph class for storing complex domains
%
% We describe a complex domain by its edges (smooth, regular
% curves) and vertices (free ends of edges or points where the ends of 
% one or more edges meets)
%
% 
    properties(SetAccess=public)
        verts
        edge2verts
        echnks
        regions
        vstruc
    end

    properties(Dependent,SetAccess=private)
        r
        d
        d2
        n
        adj
        sourceinfo
        npt
    end
    
    methods
        function obj = chunkgraph(verts,edge2verts,fchnks,cparams)
            if (nargin == 0)
                return
            end
            if (numel(verts)==0)
                return
            end
            obj.verts      = verts;
            obj.edge2verts = edge2verts;
            obj.echnks     = chunker.empty;
            
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
            obj.echnks = echnks;
            obj.vstruc = procverts(obj);
            %[regions] = findregions(obj);
            %obj.regions = regions;
            
            adjmat = edge2verts'*edge2verts;
            g = graph(adjmat);
            ccomp = conncomp(g);
            
            chnkcomp = {};
            regions = {};
            
            for i=1:max(ccomp)
                inds = find(ccomp==i);
                chnkcomp{i} = inds;
                [region_comp] = findregions(obj,inds);
                [region_comp] = findunbounded(obj,region_comp);
                regions{i} = region_comp;
            end
            
            obj.regions = regions;
        end
        function obj = set.verts(obj,val)
            classes = {'numeric'};
            validateattributes(val,classes,{})
            obj.verts = val;
        end
        function obj = set.edge2verts(obj,val)
            classes = {'numeric'};
            validateattributes(val,classes,{})
            obj.edge2verts = val;
        end    
        function obj = set.echnks(obj,val)
            classes = {'chunker'};
            validateattributes(val,classes,{})
            obj.echnks = val;
        end  
        function r = get.r(obj)
            chnk = merge(obj.echnks);
            r = chnk.r;
        end
        function d = get.d(obj)
            chnk = merge(obj.echnks);
            d = chnk.d;
        end
        function d2 = get.d2(obj)
            chnk = merge(obj.echnks);
            d2 = chnk.d2;
        end
     	function n = get.n(obj)
            chnk = merge(obj.echnks);
            n = normals(chnk);
        end
        function adj = get.adj(obj)
            chnk = merge(obj.echnks);
            adj = chnk.adj;
        end
        function npt = get.npt(obj)
            npt = 0;
            for iedge = 1:numel(obj.echnks)
            	n = size(obj.echnks(iedge).r(:,:),2);
                npt = n + npt;
            end
        end
        function sourceinfo = get.sourceinfo(obj)
            sourceinfo = [];
            ntot = obj.npt;
            rs = zeros(2,ntot);
            ds = zeros(2,ntot);
            d2s= zeros(2,ntot);
            ws = zeros(ntot,1);
            ind = 0;
            for iedge = 1:numel(obj.echnks)
                chnk = obj.echnks(iedge);
                n = chnk.npt;
                w = weights(chnk);
                ws(ind+(1:n))    = w;
                rs(:,ind+(1:n))  = chnk.r(:,:);
                ds(:,ind+(1:n))  = chnk.d(:,:);
                d2s(:,ind+(1:n)) = chnk.d2(:,:);
                ind = ind + n;
            end   
            sourceinfo.r = rs;
            sourceinfo.d = ds;
            sourceinfo.d2= d2s;
            sourceinfo.w = ws;
        end
    end   
    methods(Static)
        
    end
end
