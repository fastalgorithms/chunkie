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
% chunkgraph required properties:
%
%   verts - a (2,:) array of the coordinates of the graph vertices
%   echnks - an array of chunker objects specifying the edges of the graph
%   lrverts - a (length(echnks),2) array of the indices. lrverts(i,1) is
%      the vertex at the left end of edge i. lrverts(i,2) is the vertex at
%      the right end of edge i. A loop begins and ends at the same vertex
%
% chunkgraph derived properties:
%
%   regions - a cell array of cell arrays specifying simply connected
%      regions in the exterior of the chunkgraph. regions{j} specifies the
%      boundary of the jth region by a cell array. regions{j}{i} specifies
%      the ith connected component of the boundary of region j by an array
%      of the indices for the edges forming that component. By convention,
%      regions{1} is the unbounded connected component of the exterior.
%   vstruc - a cell array of length size(verts,2). abs(vstruc{j}) is an 
%      array indicating the edges which meet at a vertex. The sign for each
%      entry in the array vstruc{j} indicates if the edge starts (negative)
%      or ends (positive) at the vertex. If a loop meets at a vertex, the
%      edge number will appear twice in the array, with both signs.
%
%   chunkgraph merged edges properties:
%
% the properties below refer to the quantities you would get by merging all
% of the edge chunks. The points for edge 1 appear first, then edge 2, etc.
%
%   k - the order of the Legendre nodes used to discretize panels
%   nch - the total number of chunks over all edges
%   npt - the total number of points over all edges
%   r - dim x k x nch array, r(:,i,j) gives the coordinates of the ith 
%         node on the jth chunk of the chunker
%   h - nch array of scaling factors for chunks. the chunk derivatives are
%         scaled as if the coordinates in r(:,:,j) are sampled on an 
%         interval of length 2*h(j). This affects integration formulas.
%   d - dim x k x nch array, d(:,i,j) gives the time derivative of the 
%         coordinate at the ith node on the jth chunk of the chunker
%   d2 - dim x k x nch array, d(:,i,j) gives the 2nd time derivative of the 
%         coordinate at the ith node on the jth chunk of the chunker
%   n - dim x k x nch array of normals to the curve

% 
    properties(SetAccess=public)
        verts
        lrverts
        echnks
        regions
        vstruc
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
        function obj = chunkgraph(verts,lrverts,fchnks,cparams)
            if (nargin == 0)
                return
            end
            if (numel(verts)==0)
                return
            end
            obj.verts      = verts;
            obj.lrverts = lrverts;
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
            
            nvert = size(verts,2);
            if any(lrverts < 1) || any(lrverts > nvert)
                error('Incompatible vertex and connectivity (lrverts) info');
            end
            
            echnks = chunker.empty();
            for i=1:size(lrverts,1)
                i1 = lrverts(i,1);
                i2 = lrverts(i,2);
                v1 = verts(:,i1);
                v2 = verts(:,i2);
                if (numel(fchnks)<i || isempty(fchnks{i}))
                    fcurve = @(t) chnk.curves.linefunc(t,v1,v2);
                    chnkr = chunkerfunc(fcurve,cploc,pref);
                    chnkr = sort(chnkr);
                    %chnkr.vert = [v1,v2];
                    echnks(i) = chnkr;
                elseif (~isempty(fchnks{i}) && isa(fchnks{i},'function_handle'))
                    [vs,~,~] =fchnks{i}([0,1]);
                    chnkr = chunkerfunc(fchnks{i},cploc,pref);
                    chnkr = sort(chnkr);
                    v1 = verts(:,find(edge2verts(i,:)==-1));
                    v2 = verts(:,find(edge2verts(i,:)== 1));
                    r0 = vs(:,1);
                    r1 = v1;
                    scale = norm(v2-v1,'fro')/norm(vs(:,2)-vs(:,1),'fro');
                    xdfin = v2(1)-v1(1);
                    ydfin = v2(2)-v1(2);
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
    end   
    methods(Static)
        
    end
end
