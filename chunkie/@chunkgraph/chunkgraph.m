classdef chunkgraph
%CHUNKGRAPH class which describes a domain using "graph"
% of chunkers. Singular points of the geometry (or of the
% boundary value problem) are the vertices of the graph. Smooth, regular
% curves are the edges of the graph. The interiors of distinct edges should 
% not intersect. For such domains, a vertex should be introduced at the 
% point of intersection and the curves broken so that they do not intersect.
%
% chunkgraph properties:
%
%   verts - 2 x nverts array of vertex locations 
%   edgesendverts - 2 x nedges array of starting and ending vertices for
%      each edge. edgesendverts(:,i) should be NaNs for closed loops,
%      and a corresponding function handle must be provided
%   echnks - nedge x 1 array of chunker objects discretizing each edge curve
%   regions - nregions x 1 cell array specifying region information for the
%      given graph structure. A region is a connected subset of R^2 
%      specified by its bounding edges. If the jth region is simply 
%      connected then abs(regions{j}{1}) is the list of edges which 
%      comprise the boundary of region j. The sign of the edges indicate
%      the direction of traversal for that edge. The edges are ordered
%      according to an orientation based on nesting. The boundary of the
%      outermost (unbounded) region is traversed clockwise. 
%   vstruc - nverts x 1 cell array. vstruc{j}{1} gives a list of edges which
%      are incident to a vertex. vstruc{j}{2} is a list of the same length
%      consisting of +1 and -1. If +1 then the corresponding edge ends at
%      the vertex, if -1 it begins at the vertex.
%   v2emat - sparse nedge x nverts matrix, akin to a connectivity matrix.
%      the entry v2emat(i,j) is 1 if edge i ends at vertex j, it is -1 if
%      edge i starts at vertex j, and it is 2 if edge i starts and ends at
%      vertex j.
%   k - integer, number of Legendre nodes on each chunk
%   dim - integer, dimension of the ambient space in which the curve is 
%             embedded
%   In the rest, nch = \sum_{j} chunkgraph.echnks(j).nch;
%   npt - returns k*nch, the total number of points on the curve
%   r - dim x k x nch array, r(:,i,j) gives the coordinates of the ith 
%         node on the jth chunk of the chunkgraph
%   d - dim x k x nch array, d(:,i,j) gives the time derivative of the 
%         coordinate at the ith node on the jth chunk of the chunkgraph
%   d2 - dim x k x nch array, d(:,i,j) gives the 2nd time derivative of the 
%         coordinate at the ith node on the jth chunk of the chunkgraph
%   n - dim x k x nch array of normals to the curve
%   data - datadim x k x nch array of data attached to chunkgraph points 
%         this data will be refined along with the chunkgraph
%   adj - 2 x nch integer array. adj(1,j) is i.d. of the chunk that 
%         precedes chunk j in parameter space. adj(2,j) is the i.d. of the
%         chunk which follows. If adj(i,j) is 0 then that end of the chunk
%         is a free end. if adj(i,j) is negative than that end of the chunk
%         meets with other chunk ends in a vertex. 
%
% chunkgraph methods:
%   plot(obj, varargin) - plot the chunkgraph
%   quiver(obj, varargin) - quiver plot of chunkgraph points and normals
%   plot_regions(obj, iflabel) - plot the chunkgraph with region and 
%                  and edge labels
%   scatter(obj,varargin) - scatter plot of the chunkgraph nodes
%   obj = refine(obj,varargin) - refine the curve
%   wts = weights(obj) - scaled integration weights on curve
%   rn = normals(obj) - recompute normal vectors
%   obj = obj.move(r0,r1,trotat,scale) - translate, rotate, etc
%   rmin = min(obj) - minimum of coordinate values
%   rmax = max(obj) - maximum of coordinate values
%   onesmat = onesmat(obj) - matrix that corresponds to integration of a
%             density on the chunkgraph
%   rnormonesmat = normonesmat(obj) - matrix that corresponds to
%             integration of dot product of normal with vector density on
%             the chunkgraph
%   tau = tangents(obj) - unit tangents to curve
%   kappa = signed_curvature(obj) - get signed curvature along curve
%   obj = makedatarows(obj,nrows) - add nrows rows to the data storage.
%   rflag = datares(obj,opts) - check if data in data rows is resolved
%
%   To add:
%     flagnear
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
% 
    properties(SetAccess = public)
        verts
        edgesendverts
        echnks
        regions
        vstruc
        v2emat
   
        r
        d
        d2
        n
        wts
        adj
        sourceinfo

        data
    end

    properties(Dependent, SetAccess=private)
        npt
        k
        dim
        datadim
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
            obj.v2emat        = build_v2emat(obj);
            obj.echnks        = chunker.empty;

            if nargin < 3
                fchnks = [];
            end
            
            if isa(fchnks,"function_handle")
                fchnks0 = fchnks;
                fchnks = cell(nedge,1);
                for j = 1:nedge
                    fchnks{j} = fchnks0;
                end
            end

            if (nargin < 4)
                cparams = [];
            end
            msg = "CHUNKGRAPH: cparams must be struct or nedge cell array of structs";
            assert(isempty(cparams) || isstruct(cparams) || (iscell(cparams) && length(cparams) == nedge),msg);

            if nargin < 5
                pref = chunkerpref();
            else
                pref = chunkerpref(pref);
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

                    if (isnan(i1) || isnan(i2))
                        msg = "CHUNKGRAPH:CONSTRUCTOR: cannot create a " + ...
                              "smooth curve without a function handle." + ...
                              " fchnks{iedge} must be provided if either" + ...
                              "vertex of edge end is NaN";
                        error(msg);
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
                    i1 = edgesendverts(1,i);
                    i2 = edgesendverts(2,i);

                    if isnan(i1) || isnan(i2)
                        cploc.ifclosed = true;
                    end

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
                    if isfield(cploc,'maxchunklen')
                        if ~isnan(i1) && ~isnan(i2)
                            if i1 ~= i2 
                                vfin0 = verts(:,i1);
                                vfin1 = verts(:,i2);
                                scale = norm(vfin1-vfin0,'fro')/norm(vs(:,2)-vs(:,1),'fro');
                                cploc.maxchunklen = cploc.maxchunklen/scale;
                            end
                        end
                    end
                    
                    chnkr = chunkerfunc(fchnks{i}, cploc, pref);
                    chnkr = sort(chnkr);
                    if ~isnan(i1) && ~isnan(i2)
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
                    end

                    echnks(i) = chnkr;
                end
            end
            obj.echnks = echnks;
            obj.vstruc = procverts(obj);
            obj.wts = weights(obj);
           
            
            obj.regions = findregions(obj);

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
        function k = get.k(obj)
            k = size(obj.r,2);
        end
        function dim = get.dim(obj)
            dim = size(obj.r,1);
        end
        function datadim = get.datadim(obj)
            datadim = size(obj.data,1);
        end

        function sourceinfo = get.sourceinfo(obj)
            sourceinfo = [];
            chnkrtotal = merge(obj.echnks);
            
            sourceinfo.r = chnkrtotal.r(:,:);
            sourceinfo.n = chnkrtotal.n(:,:);
            sourceinfo.d = chnkrtotal.d(:,:);
            sourceinfo.d2 = chnkrtotal.d2(:,:);
            sourceinfo.w = chnkrtotal.wts(:);

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
        obj = refine(obj, opts)
        obj = balance(obj)
        obj = move(obj, r0, r1, trotat, scale)
        rmin = min(obj)
        rmax = max(obj)
        plot(obj, varargin)
        quiver(obj, varargin)
        plot_regions(obj, iflabel)
        wts = weights(obj)
        rnorm = normals(obj)
        onesmat = onesmat(obj)
        rnormonesmat = normonesmat(obj)
        tau = tangents(obj)
        kappa = signed_curvature(obj)
        obj = makedatarows(obj, nrows)
        scatter(obj, varargin)
        rres = datares(obj, opts)

        
    end

    methods(Static)
    end
end
