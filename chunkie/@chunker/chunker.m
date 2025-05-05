classdef chunker
%CHUNKER class which describes a curve divided into chunks (or "panels"). 
%
% On each chunk the curve is represented by the values of its position, 
% first and second derivatives by scaled Legendre nodes.
%
% chunker properties:
%   k - integer, number of Legendre nodes on each chunk
%   nch - integer, number of chunks that make up the curve
%   dim - integer, dimension of the ambient space in which the curve is 
%             embedded
%   npt - returns k*nch, the total number of points on the curve
%   r - dim x k x nch array, r(:,i,j) gives the coordinates of the ith 
%         node on the jth chunk of the chunker
%   d - dim x k x nch array, d(:,i,j) gives the time derivative of the 
%         coordinate at the ith node on the jth chunk of the chunker
%   d2 - dim x k x nch array, d(:,i,j) gives the 2nd time derivative of the 
%         coordinate at the ith node on the jth chunk of the chunker
%   n - dim x k x nch array of normals to the curve
%   data - datadim x k x nch array of data attached to chunker points 
%         this data will be refined along with the chunker
%   adj - 2 x nch integer array. adj(1,j) is i.d. of the chunk that 
%         precedes chunk j in parameter space. adj(2,j) is the i.d. of the
%         chunk which follows. if adj(i,j) is 0 then that end of the chunk
%         is a free end. if adj(i,j) is negative than that end of the chunk
%         meets with other chunk ends in a vertex. the specific negative
%         number acts as an i.d. of that vertex.
%
% chunker methods: 
%   chunker(p,t,w) - construct an empty chunker with given preferences and
%       precomputed Legendre nodes/weights (optional)
%   obj = addchunk(obj,nchadd) - add nchadd chunks to the structure 
%       (initialized with zeros)
%   obj = obj.move(r0,r1,trotat,scale) - translate, rotate, etc
%   obj = makedatarows(obj,nrows) - add nrows rows to the data storage.
%   [obj,info] = sort(obj) - sort the chunks so that adjacent chunks are
%        stored sequentially
%   [rn,dn,d2n,dist,tn,ichn] = nearest(obj,ref,ich,opts,u,xover,aover) -
%        find nearest point on chunker to ref
%   obj = reverse(obj) - reverse chunk orientation
%   rmin = min(obj) - minimum of coordinate values
%   rmax = max(obj) - maximum of coordinate values
%   wts = weights(obj) - scaled integration weights on curve
%   obj.n = normals(obj) - recompute normal vectors
%   onesmat = onesmat(obj) - matrix that corresponds to integration of a
%             density on the curve
%   rnormonesmat = normonesmat(obj) - matrix that corresponds to
%             integration of dot product of normal with vector density on
%             the curve
%   plot(obj,varargin) - plot the chunker curve
%   plot3(obj,idata,varargin) - 3D plot of the curve and one row of the
%             data storage
%   quiver(obj,varargin) - quiver plot of the chnkr points and normals
%   scatter(obj,varargin) - scatter plot of the chnkr nodes
%   tau = taus(obj) - unit tangents to curve
%   obj = refine(obj,varargin) - refine the curve
%   a = area(obj) - for a closed curve, area inside 
%   s = arclength(obj) - get values of arclength along curve
%   kappa = signed_curvature(obj) - get signed curvature along curve
%   rflag = datares(obj,opts) - check if data in data rows is resolved
%   [rc,dc,d2c] = exps(obj) - get expansion coefficients for r,d,d2
%   ier = checkadjinfo(obj) - checks that the adjacency info of the curve
%              is consistent
%   [inds,adjs,info] = sortinfo(obj) - attempts to sort the curve and finds
%              number of connected components, etc
%   [re,taue] = chunkends(obj,ich) - get the endpoints of chunks
%   flag = flagnear(obj,pts,opts) - flag points near the boundary


% author: Travis Askham (askhamwhat@gmail.com)

    properties(Access=private)
        verttol
    end
    properties(Dependent,Access=public)
        r
        d
        d2
        adj
        n
        wts
        data
    end
    properties(Hidden, Access=public)
        rstor
        dstor
        d2stor
        adjstor
        nstor
        wtsstor
        datastor
        hasdata
    end
    properties(SetAccess=private)
        nch
    end
    properties(Hidden, SetAccess=private)
        nchmax
        nchstor
        vert
        wstor
        tstor
    end
    properties(Dependent,SetAccess=private)
        k
        dim
        npt
        datadim
    end
    properties(Hidden,Dependent,SetAccess=private)
        nvert
        vertdeg
    end
    
    methods
        function obj = chunker(p,t,w)
        %CHUNKER construct an empty chunker with some default settings
        %
        % syntax: chnkr = chunker(p,t,w);
        %
        % optional input: 
        %   p - struct of preferences 
        %        p.k - integer, order to be used on chunks (16)
        %        p.dim - integer, dimension of ambient space (2)
        %   t, w - arrays of Legendre nodes and weights of order p.k
        % 
            if nargin < 1 || isempty(p)
                p = chunkerpref();
            else
                p = chunkerpref(p);
            end
            k = p.k;
            assert(k >= 2,'CHUNKER: order k of panels must be at least 2');
            if nargin < 3
                [obj.tstor,obj.wstor] = lege.exps(k);
            else
                assert(length(t)==k && length(w)==k,...
                    'CHUNKER: precomputed Legendre nodes appear to be wrong order');
                obj.tstor = t;
                obj.wstor = w;
            end
            dim = p.dim;
            nchmax = p.nchmax;
            nchstor = p.nchstor;
            obj.verttol = p.verttol;
            obj.nchmax = nchmax;
            obj.nchstor = nchstor;
            obj.nch = 0;
            obj.rstor = zeros(dim,k,nchstor);
            obj.dstor = zeros(dim,k,nchstor);
            obj.nstor = zeros(dim,k,nchstor);
            obj.wtsstor = zeros(k,nchstor);
            obj.d2stor = zeros(dim,k,nchstor);
            obj.adjstor = zeros(2,nchstor);
            obj.vert = {};
            obj.hasdata = false;
            obj.datastor = [];
            
        end
        
        function r = get.r(obj)
            r = obj.rstor(:,:,1:obj.nch);
        end
        function d = get.d(obj)
            d = obj.dstor(:,:,1:obj.nch);
        end
        function d2 = get.d2(obj)
            d2 = obj.d2stor(:,:,1:obj.nch);
        end
        function adj = get.adj(obj)
            adj = obj.adjstor(:,1:obj.nch);
        end
        function n = get.n(obj)
            n = obj.nstor(:,:,1:obj.nch);
        end
        function wts = get.wts(obj)
            wts = obj.wtsstor(:,1:obj.nch);
        end
        function data = get.data(obj)
            data = obj.datastor(:,:,1:(obj.nch*obj.hasdata));
        end
        function obj = set.data(obj,val)
            obj.datastor(:,:,1:(obj.nch*obj.hasdata)) = val;
        end
        function obj = set.r(obj,val)
            obj.rstor(:,:,1:obj.nch)=val;
        end
        function obj = set.d(obj,val)
            obj.dstor(:,:,1:obj.nch) = val;
        end
        function obj = set.wts(obj,val)
            obj.wtsstor(:,1:obj.nch) = val;
        end        
        function obj = set.n(obj,val)
            obj.nstor(:,:,1:obj.nch) = val;
        end        
        function obj = set.d2(obj,val)
            obj.d2stor(:,:,1:obj.nch) = val;
        end
        function obj = set.adj(obj,val)
            obj.adjstor(:,1:obj.nch) = val;
        end
        function k = get.k(obj)
            k = size(obj.rstor,2);
        end
        function dim = get.dim(obj)
            dim = size(obj.rstor,1);
        end
        function nvert = get.nvert(obj)
            nvert = numel(obj.vert);
        end
        function datadim = get.datadim(obj)
            datadim = size(obj.data,1);
        end
        function npt = get.npt(obj)
            npt = obj.k*obj.nch;
        end
        function vertdeg = get.vertdeg(obj)
            vertdeg = cellfun(@numel,obj.vert);
        end
            
        function obj = addchunk(obj,nchadd)
            if nargin < 2
                nchadd = 1;
            end
            assert(and(isnumeric(nchadd),nchadd > 0), ...
                'nchadd must be positive integer');
            assert(obj.nch+nchadd <= obj.nchmax, ...
                'adding chunks would exceed maximum chunker length');
            while obj.nch + nchadd > obj.nchstor
                % double size until sufficient storage available
                obj = resize(obj,max(min(2*obj.nchstor,obj.nchmax),...
                    min(obj.nch+nchadd,obj.nchmax)));
            end
            obj.nch = obj.nch+nchadd;
            assert(obj.nch <= min(obj.nchstor,obj.nchmax),'something went wrong');
        end
        
        function obj = resize(obj,nchstornew)
            assert(nchstornew >= obj.nch, ...
                ['error new length is less than number of chunks\n', ...
                'info would be lost. coarsen/delete chunks first to shrink']);
            assert(nchstornew <= obj.nchmax, ...
                ['new storage exceeds maximum storage.\n', ...
                'perhaps something went wrong. otherwise, increase nchmax']);
            rtemp = obj.r;
            dtemp = obj.d;
            ntemp = obj.n;
            wtstemp = obj.wts;
            d2temp = obj.d2;
            adjtemp = obj.adj;
            datatemp = obj.data;
            obj.rstor = zeros(obj.dim,obj.k,nchstornew);
            obj.dstor = zeros(obj.dim,obj.k,nchstornew);
            obj.nstor = zeros(obj.dim,obj.k,nchstornew);
            obj.wtsstor = zeros(obj.k,nchstornew);
            obj.d2stor = zeros(obj.dim,obj.k,nchstornew);
            obj.adjstor = zeros(2,nchstornew);
            obj.datastor = zeros(obj.datadim,obj.k,nchstornew);
            obj.r = rtemp;
            obj.d = dtemp;
            obj.n = ntemp;
            obj.wts = wtstemp;
            obj.d2 = d2temp;
            obj.adj = adjtemp;
            obj.data = datatemp;            
            obj.nchstor = nchstornew;
        end
        
        function obj = makedatarows(obj,nrows)
            if (nrows > 0)
                datatemp = obj.data;
                datadimold = obj.datadim;
                obj.datastor = zeros(datadimold+nrows,obj.k,obj.nchstor);
                obj.data(1:datadimold,:) = datatemp(:,:);
                obj.hasdata = true;
            else
                if (nrows < 0)
                    warning('attempted to add negative rows, doing nothing');
                end
            end
        end
        
        function obj = clearverts(obj)
            obj.vert = {};
        end
        function obj = addvert(obj,newvert,toleft)
        %ADDVERT add vertex label to chunker object 
        %
            if (numel(newvert) == 0); warning('chunker:emptyvert',...
                    'no vertex specified, doing nothing'); end
            if (numel(newvert) == 1); warning('chunker:badvert',...
              'free ends are not vertices in chunkers, doing nothing'); end
            if (numel(newvert) > 1)
                assert(all(and(1 <= newvert,newvert <= obj.nch)),...
                    'vertex indices must correspond to existing chunks');
                ends = chunkends(obj,newvert(:));
                if (nargin < 3)
                    lends = ends(:,1,:); rends = ends(:,2,:);
                    lperm = permute(lends,[1,3,2]); 
                    rperm = permute(rends,[1,3,2]);
                
                    lldist = squeeze(sqrt(sum( (lends - lperm).^2, 1)));
                    rldist = squeeze(sqrt(sum( (lends - rperm).^2, 1)));
                    rrdist = squeeze(sqrt(sum( (rends - rperm).^2, 1)));
                
                    ldist = min(lldist,rldist.'); ldist = sum(ldist,2);
                    rdist = min(rldist,rrdist); rdist = sum(rdist,2);
                    toleft = ldist < rdist;
                end
                
                lens = chunklen(obj,newvert); 
                maxlen = max(lens(:));
                tol = obj.verttol; tol = max(tol,maxlen*tol);
                
                vertends = zeros(obj.dim,numel(newvert));
                vertends(:,toleft) = ends(:,1,toleft);
                vertends(:,~toleft) = ends(:,2,~toleft);
                
                vertends = vertends - mean(vertends,2);
                errvertends = norm(vertends,'fro');
                if (errvertends > tol)
                    warning('chunkie:badvertchunks', ...
                     [sprintf('chunkends far from centroid %5.2e\n',...
                     errvertends),...
                     'consider resetting verttol for chunker']);
                end
                
                nvert1 = obj.nvert;
                nvert1 = nvert1+1;
                obj.adjstor(1,newvert(toleft)) = -nvert1;
                obj.adjstor(2,newvert(~toleft)) = -nvert1;
                obj.vert{nvert1} = newvert;
                
            end
        end    
            
        
        function obj = cleardata(obj)
            obj.hasdata = false;
            obj.datastor = [];
        end
        
        [obj2,f2] = upsample(obj,kup,f)
        [obj,info] = sort(obj)
        [rn,dn,d2n,dist,tn,ichn] = nearest(obj,ref,ich,opts,u,xover,aover)
        obj = reverse(obj)
        obj = move(obj,r0,r1,trotat,scale)
        rmin = min(obj)
        rmax = max(obj)
        whts = whts(obj)
        wts = weights(obj)
        rnorm = normals(obj)
        onesmat = onesmat(obj)
        rnormonesmat = normonesmat(obj)
        df = diff(obj,f,ndim)
        plot(obj,varargin)
        plot3(obj,idata,varargin)
        quiver(obj,varargin)
        scatter(obj,varargin)
        tau = taus(obj)
        obj = refine(obj,varargin)
        a = area(obj)
        s = arclength(obj)
        rflag = datares(obj,opts)
        [rc,dc,d2c] = exps(obj)
        ier = checkadjinfo(obj)
        [inds,adjs,info] = sortinfo(obj)
        [re,taue] = chunkends(obj,ich)
        flag = flagnear(obj,pts,opts)
        kappa = signed_curvature(obj)
        obj = plus(v,obj)
        obj = mtimes(A,obj)
    end
    methods(Static)
        obj = chunkerfunc(fcurve,varargin)
        obj = chunkerpoly(verts,varargin)
        function lvlrfac = lvlrfacdefault()
            lvlrfac = 2.1;
        end
    end
end
