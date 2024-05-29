classdef chunkgraph
%CHUNKGRAPH chunk graph class for storing complex domains
%
% We describe a complex domain by a planar graph. Smooth boundary components
% (discretized by chunkers) specify the edges of the graph and the
% coordinates of free ends of the edges or points where the ends of multiple
% edges meet specify the vertices. The interiors of distinct edges are not 
% allowed to intersect (i.e. intersecting curves should be split at the
% point where they intersect. 
%
% The chunkgraph is specified by an array of chunkers called echnks (the
% edges), an array of points called verts (the vertices), and an array of
% indices specifying the vertices at the left and right ends of any 
%
% 
    properties(SetAccess=public)
        verts
        edgesendverts
        echnks
        regions
        vstruc
        v2emat
    end

    properties(SetAccess=public)
        r
        d
        d2
        n
        wts
        adj
        sourceinfo
        npt
    end
    
    methods
        function obj = chunkgraph(verts, edgesendverts, fchnks, cparams, pref)
            if (nargin == 0)
                return
            end
            if (numel(verts)==0)
                return
            end
            obj.verts = verts;

            nverts = size(verts(:,:),2);
            assert(nverts == size(edgesendverts,2),'edge specification not compatible with number of vertices');
            if (size(edgesendverts,1) ~= 2)
                nedge = size(edgesendverts,1);
                edgevertends_new = zeros(2,nedge);
                nedge = size(edgesendverts,1);
                for i = 1:nedge
                    edgevertends_new(1,i) = find(edgesendverts(i,:) == -1);
                    edgevertends_new(2,i) = find(edgesendverts(i,:) == 1);
                end
                edgesendverts = edgevertends_new;
                assert(all(edgesendverts(:) ~= 0),'edge specification had an error');
            end 

            obj.edgesendverts = edgesendverts;
            obj.v2emat = build_v2emat(obj);
            obj.echnks     = chunker.empty;

            if nargin < 3
                fchnks = [];
            end
            
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
            
            if nargin <= 4
                pref = [];
                pref.nchmax = 10000;
                pref.k = 16;
            end
            
            if ~isfield(pref, 'nchmax')
                pref.nchmax = 10000;
            end
            
            
            echnks = chunker.empty();
            for i=1:size(edgesendverts,2)
                if (numel(fchnks)<i || isempty(fchnks{i}))
                    i1 = edgesendverts(1,i);
                    i2 = edgesendverts(2,i);
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
                    i1 = edgesendverts(1,i);
                    i2 = edgesendverts(2,i);
                    vfin0 = verts(:,i1);
                    vfin1 = verts(:,i2);
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
            obj.wts = weights(obj);
            %[regions] = findregions(obj);
            %obj.regions = regions;
           
            g = graph(edgesendverts(1,:),edgesendverts(2,:));
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
            
            gmat = zeros(numel(regions));
            
            for ii=1:numel(regions)
               rgna = regions{ii};
               ilist = [];
               for jj=1:numel(regions)
                   if (ii ~=jj)
                        rgnb = regions{jj};
                        [inc] = regioninside(obj,rgnb,rgna);
                        if (inc)
                            ilist = [ilist,jj];
                        end
                   end
                   gmat(ii,ilist) = 1;
                   gmat(ilist,ii) = 1;
               end
               imin = min(ilist);   
            end    
            
            ccomp_reg = conncomp(graph(gmat));
            [s,inds] = sort(ccomp_reg);
            regions = regions(inds);
            
            for ii = 1:numel(s)
                si = s(ii);
                for jj=1:(numel(s)-1)
                    sj = s(jj);
                    if (si == sj)
                        rgna = regions{jj};
                        rgnb = regions{jj+1};
                        [inc] = regioninside(obj,rgna,rgnb);
                        if (inc)
                            regions([jj,jj+1])= regions([jj+1,jj]);
                        end
                    end
                end    
            end

            
            rgns = regions;
            rgnso= {};
            
            for ii=1:max(s)
                inds = find(s==ii);
                rgnout = rgns{inds(1)};
                for jj=2:numel(inds)
                    indj = inds(jj);
                    [rgnout] = mergeregions(obj,rgnout,rgns{indj});
                end
                rgnso{ii} = rgnout;
            end    
            
            regions = rgnso;
            rgns = regions;
            rgnout = rgns{1};
            if (numel(rgns)>1)
                rgn2 = rgns{2};
                [rgnout] = mergeregions(obj,rgnout,rgn2);
                for ii=3:numel(rgns)
                    rgn2 = rgns{ii};
                    [rgnout] = mergeregions(obj,rgnout,rgn2);
                end
            end
            
            regions = rgnout;
            
            obj.regions = regions;
        end
        function obj = set.verts(obj,val)
            classes = {'numeric'};
            validateattributes(val,classes,{})
            obj.verts = val;
        end
        function obj = set.edgesendverts(obj,val)
            classes = {'numeric'};
            validateattributes(val,classes,{})
            obj.edgesendverts = val;
        end    
        function obj = set.echnks(obj,val)
            classes = {'chunker'};
            validateattributes(val,classes,{})
            obj.echnks = val;
        end  
        function obj = set.wts(obj,val)
            obj.wts = val;
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
            n = chnk.n;
        end
        function wts = get.wts(obj)
            chnk = merge(obj.echnks);
            wts = chnk.wts;
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
                w = chnk.wts;
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
        function inds = edgeinds(obj,edgelist)
            ladr = cumsum([1,obj.echnks.npt]);
            inds = [];
            for j = 1:length(edgelist)
                ej = edgelist(j);
                inds = [inds,ladr(ej):(ladr(ej+1)-1)];
            end
        end

        % defined in other files 
        spmat = build_v2emat(obj)
    end

    methods(Static)
    end
end
