classdef chunker
%CHUNKER class which describes a curve divided 
% up into chunks. On each chunk the curve is represented by 
% the values of its position, first and second derivatives in 
% parameter space on a Legendre grid
    properties(Access=private)
        rstor
        dstor
        d2stor
        adjstor
        hstor
        datastor
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
    end
    properties(Dependent,SetAccess=private)
        k
        dim
        npt
        datadim
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
            obj.nchmax = nchmax;
            obj.nchstor = nchstor;
            obj.nch = 0;
            obj.rstor = zeros(dim,k,nchstor);
            obj.dstor = zeros(dim,k,nchstor);
            obj.d2stor = zeros(dim,k,nchstor);
            obj.adjstor = zeros(2,nchstor);
            obj.hstor = zeros(nchstor,1);
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
        function h = get.h(obj)
            h = obj.hstor(1:obj.nch);
        end
        function data = get.data(obj)
            data = obj.datastor(:,1:(obj.nch*obj.hasdata));
        end
        function obj = set.data(obj,val)
            obj.datastor(:,1:(obj.nch*obj.hasdata)) = val;
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
            k = size(obj.r,2);
        end
        function dim = get.dim(obj)
            dim = size(obj.r,1);
        end
        function datadim = get.datadim(obj)
            datadim = size(obj.data,1);
        end
        function npt = get.npt(obj)
            npt = obj.k*obj.nch;
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
            obj.rstor = zeros(obj.dim,obj.k,nchstornew);
            obj.dstor = zeros(obj.dim,obj.k,nchstornew);
            obj.d2stor = zeros(obj.dim,obj.k,nchstornew);
            obj.adjstor = zeros(2,nchstornew);
            obj.hstor = zeros(nchstornew,1);
            obj.r = rtemp;
            obj.d = dtemp;
            obj.d2 = d2temp;
            obj.adj = adjtemp;
            obj.h = htemp;
            obj.nchstor = nchstornew;
        end
        
        function obj = makedatarows(obj,nrows)
            assert(nrows > 0, ...
                'error must add a positive number of rows to data');
            datatemp = obj.data;
            datadimold = obj.datadim;
            obj.datastor = zeros(datadimold+nrows,obj.k,obj.nchstor);
            obj.datastor(1:datadimold,:) = datatemp(:,:);
        end
        
        [obj,ifclosed] = sort(obj)
        obj = reverse(obj)
        rmin = min(obj)
        rmax = max(obj)
        whts = whts(obj)
        rnorm = normals(obj)
        onesmat = onesmat(obj)
        rnormonesmat = normonesmat(obj)
        df = diff(obj,f,ndim)
        plot(obj,varargin)
        quiver(obj,varargin)
        scatter(obj,varargin)
        tau = taus(obj)
        obj = refine(obj,varargin)
        a = area(obj)
    end
    methods(Static)
        obj = chunkfunc(fcurve,varargin)
    end
end