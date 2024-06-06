classdef chunkgraph
%CHUNKGRAPH chunk graph class which describes a domain using "graph"
% terminology. Singular points of the geometry (or of the
% boundary value problem) are the vertices of the graph. Smooth, regular
% curves are the edges of the graph. A vertex must be the end point of at
% least one edge. The interiors of distinct edges should not intersect. For
% such domains, a vertex should be introduced at the point of intersection
% and the curves broken so that they do not intersect.
%
% Syntax:
%
%      cg = chunkgraph(verts,edgesendverts,fchnks,cparams,pref)
%
% will construct a chunkgraph.
%
% Input:
%   verts - (2 x nverts) array specifying vertex locations, verts(:,j) is
%     the jth vertex.
%   edgesendverts - (2 x nedges) array specifying starting and ending
%     vertices of an edge. 
% Optional input:
%   fchnks - (nedges x 1) cell array of function handles, specifying a
%     smooth curve to connect the given vertices. fchnk{j} should be a
%     function in the format expected by CHUNKERFUNC but with the default
%     that fchnk{j} is a function from [0,1] to the curve. If the specified
%     curve does not actually connect the given vertices, the curve will be
%     translated, rotated, and scaled to connect them. If the specified
%     curve should be a loop, it is only translated to start at the correct
%     vertex. 
%   cparams - struct or (nedges x 1) cell array of structs specifying curve
%     parameters in the format expected by CHUNKERFUNC. 
%   pref - struct specifying CHUNKER preferences. 
%  
%   Note:
%
%   if these inputs are specified, then 
%               echnks(j) = chunkerfun(fchnks{j},cparams{j},pref)
%   if no function handle is specified, then echnks(j) will be the
%   straight line connecting the given vertices. 
%
%  Class properties:
%
%  verts = 2 x nverts array of vertex locations 
%  edgesendverts = 2 x nedges array of starting and ending vertices for
%      each edge 
%  echnks = nedge x 1 array of chunker objects discretizing each edge curve
%  regions = nregions x 1 cell array specifying region information for the
%      given graph structure. A region is a connected subset of R^2 
%      specified by its bounding edges. If the jth region is simply 
%      connected then abs(regions{j}{1}) is the list of edges which 
%      comprise the boundary of region j. The sign of the edges indicate
%      the direction of traversal for that edge. The edges are ordered
%      according to an orientation based on nesting. The boundary of the
%      outermost (unbounded) region is traversed clockwise. 
%  vstruc = nverts x 1 cell array. vstruc{j}{1} gives a list of edges which
%      are incident to a vertex. vstruc{j}{2} is a list of the same length
%      consisting of +1 and -1. If +1 then the corresponding edge ends at
%      the vertex, if -1 it begins at the vertex.
%  v2emat = sparse nedge x nverts matrix, akin to a connectivity matrix.
%      the entry v2emat(i,j) is 1 if edge i ends at vertex j, it is -1 if
%      edge i starts at vertex j, and it is 2 if edge i starts and ends at
%      vertex j.
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
            %CHUNKGRAPH constructor. Documented above.
            %
            if (nargin == 0)
                return
            end
            if (numel(verts)==0)
                return
            end
            obj.verts = verts;

            nverts = size(verts(:,:),2);
            if (size(edgesendverts,1) ~= 2 || any(edgesendverts(:) <= 0))
                assert(nverts == size(edgesendverts,2),'edge specification not compatible with number of vertices');

                nedge = size(edgesendverts,1);
                edgevertends_new = zeros(2,nedge);
                nedge = size(edgesendverts,1);
                for i = 1:nedge
                    edgevertends_new(1,i) = find(edgesendverts(i,:) == -1);
                    edgevertends_new(2,i) = find(edgesendverts(i,:) == 1);
                end
                edgesendverts = edgevertends_new;
                assert(all(edgesendverts(:) ~= 0),'edge specification had an error');
            else
                nedge = size(edgesendverts,2);
            end 

            obj.edgesendverts = edgesendverts;
            obj.v2emat = build_v2emat(obj);
            obj.echnks     = chunker.empty;

            if nargin < 3
                fchnks = [];
            end
            
            if (nargin < 4)
                cparams = [];
            end
            msg = "CHUNKGRAPH: cparams must be struct or nedge cell array of structs";
            assert(isempty(cparams) || isstruct(cparams) || (iscell(cparams) && length(cparams) == nedge),msg);

            if nargin < 5
                pref = [];
            end
            
            echnks = chunker.empty();
            for i=1:size(edgesendverts,2)
                if (numel(fchnks)<i || isempty(fchnks{i}))
                    if iscell(cparams)
                        cploc = cparams{i};
                    else
                        cploc = cparams;
                    end
                    % set cploc.ifclosed in a way that makes sense
                    cploc.ifclosed = false;
                    % chunkgraph edges need at least 4 chunks
                    if isfield(cploc,'nchmin')
                        cploc.nchmin = max(4,cploc.nchmin);
                    else
                        cploc.nchmin = 4;
                    end
                    % line function uses these ta and tb
                    cploc.ta = 0;
                    cploc.tb = 1;
                    
                    i1 = edgesendverts(1,i);
                    i2 = edgesendverts(2,i);
                    if (i1 == i2)
                        msg = "chunkgraph constructor: connecting a " + ...
                            "vertex to itself by a line " + ...
                            "may have unexpected behavior";
                        warning(msg);
                    end
                    v1 = verts(:,i1);
                    v2 = verts(:,i2);
                    fcurve = @(t) chnk.curves.linefunc(t,v1,v2);
                    chnkr = chunkerfunc(fcurve,cploc,pref);
                    chnkr = sort(chnkr);
                    %chnkr.vert = [v1,v2];
                    echnks(i) = chnkr;
                elseif (~isempty(fchnks{i}) && isa(fchnks{i},'function_handle'))
                    if iscell(cparams)
                        cploc = cparams{i};
                    else
                        cploc = cparams;
                    end
                    % set cploc.ifclosed in a way that makes sense
                    cploc.ifclosed = false;
                    % chunkgraph edges need at least 4 chunks
                    if isfield(cploc,'nchmin')
                        cploc.nchmin = max(4,cploc.nchmin);
                    else
                        cploc.nchmin = 4;
                    end
                    
                    ta = 0; tb = 1;
                    if isfield(cploc,'ta')
                        ta = cploc.ta; 
                    else
                        cploc.ta = ta;
                    end
                    if isfield(cploc,'tb')
                        tb = cploc.tb; 
                    else
                        cploc.tb = tb;
                    end
                    vs =fchnks{i}([ta,tb]);
                    chnkr = chunkerfunc(fchnks{i},cploc,pref);
                    chnkr = sort(chnkr);
                    i1 = edgesendverts(1,i);
                    i2 = edgesendverts(2,i);
                    if i1 ~= i2
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
                    else
                        vfin = verts(:,i1);
                        if norm(vs(:,1)-vs(:,2))/max(abs(chnkr.r(:))) > 1e-14
                            msg = "chunkgraph constructor: edge " + ...
                            num2str(i) + " is defined as a loop but " + ...
                            "curve defining the edge does not reconnect " + ...
                            "to high precision";
                            warning(msg);
                        end
                        chnkr = move(chnkr,zeros(size(vs(:,1))),vfin-vs(:,1),0,1);
                    end

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

            obj = balance(obj);
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
        obj = refine(obj,opts)
        obj = balance(obj)
    end

    methods(Static)
    end
end
