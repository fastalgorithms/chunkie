classdef chunker
%CHUNKER class which describes a curve divided into chunks (panels). 
% On each chunk the curve is represented by the values of its position, 
% first and second derivatives in parameter space on a Legendre grid.
%

% author: Travis Askham (askhamwhat@gmail.com)

    properties(Access=private)
        rstor
        dstor
        d2stor
        adjstor
        hstor
        datastor
        verttol
        wstor
        tstor
    end
    properties(Dependent,Access=public)
        r
        d
        d2
        adj
        h
        data
    end
    properties(SetAccess=private)
        nchmax
        nch
        nchstor
        hasdata
        vert
    end
    properties(Dependent,SetAccess=private)
        k
        dim
        npt
        datadim
        nvert
        vertdeg
    end
    
    methods
        function obj = chunker(p)
            if nargin < 1
                p = chunkerpref();
            else
                p = chunkerpref(p);
            end
            k = p.k;
            dim = p.dim;
            nchmax = p.nchmax;
            nchstor = p.nchstor;
            obj.verttol = p.verttol;
            obj.nchmax = nchmax;
            obj.nchstor = nchstor;
            obj.nch = 0;
            obj.rstor = zeros(dim,k,nchstor);
            obj.dstor = zeros(dim,k,nchstor);
            obj.d2stor = zeros(dim,k,nchstor);
            obj.adjstor = zeros(2,nchstor);
            obj.hstor = zeros(nchstor,1);
            obj.vert = {};
            obj.hasdata = false;
            obj.datastor = [];
            [obj.tstor,obj.wstor] = lege.exps(k);
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
        function h = get.h(obj)
            h = obj.hstor(1:obj.nch);
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
        function obj = set.d2(obj,val)
            obj.d2stor(:,:,1:obj.nch) = val;
        end
        function obj = set.adj(obj,val)
            obj.adjstor(:,1:obj.nch) = val;
        end
        function obj = set.h(obj,val)
            obj.hstor(1:obj.nch) = val;
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
            d2temp = obj.d2;
            adjtemp = obj.adj;
            htemp = obj.h;
            datatemp = obj.data;
            obj.rstor = zeros(obj.dim,obj.k,nchstornew);
            obj.dstor = zeros(obj.dim,obj.k,nchstornew);
            obj.d2stor = zeros(obj.dim,obj.k,nchstornew);
            obj.adjstor = zeros(2,nchstornew);
            obj.hstor = zeros(nchstornew,1);
            obj.datastor = zeros(obj.datadim,obj.k,nchstornew);
            obj.r = rtemp;
            obj.d = dtemp;
            obj.d2 = d2temp;
            obj.adj = adjtemp;
            obj.h = htemp;
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
                
                lens = chunklen(obj,newvert,obj.wstor); 
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
                     errvertends),
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
        
        [obj,info] = sort(obj)
        [rn,dn,d2n,dist,tn,ichn] = nearest(obj,ref,ich,opts,u,xover,aover)
        obj = reverse(obj)
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
    end
    methods(Static)
        obj = chunkerfunc(fcurve,varargin)
    end
end
