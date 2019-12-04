classdef trapper
%TRAPPER class which describes a curve based on
% equispaced sampling in some parameter space

    properties(Access=private)
        rstor
        dstor
        d2stor
        datastor
    end
    properties(Dependent,Access=public)
        r
        d
        d2
        data
    end
    properties(Access=public)
        h 
    end
    properties(SetAccess=private)
        npt
        nptstor
        nptmax
        hasdata
    end
    properties(Dependent,SetAccess=private)
        dim
        datadim
    end
    
    methods
        function obj = trapper(p)
            if nargin < 1
                p = trapperpref();
            else
                p = trapperpref(p);
            end
            dim = p.dim;
            nptstor = p.nptstor;
            nptmax = p.nptmax;
            obj.rstor = zeros(dim,nptstor);
            obj.dstor = zeros(dim,nptstor);
            obj.d2stor = zeros(dim,nptstor);
            obj.npt = 0;
            obj.nptmax = nptmax;
            obj.hasdata = false;
            obj.datastor = [];
        end
        
        function r = get.r(obj)
            r = obj.rstor(:,1:obj.npt);
        end
        function d = get.d(obj)
            d = obj.dstor(:,1:obj.npt);
        end
        function d2 = get.d2(obj)
            d2 = obj.d2stor(:,1:obj.npt);
        end
        function h = get.h(obj)
            h = obj.h;
        end
        function data = get.data(obj)
            data = obj.datastor(:,1:(obj.npt*obj.hasdata));
        end
        function obj = set.data(obj,val)
            obj.datastor(:,1:(obj.npt*obj.hasdata)) = val;
        end
        function obj = set.r(obj,val)
            obj.rstor(:,1:obj.npt)=val;
        end
        function obj = set.d(obj,val)
            obj.dstor(:,1:obj.npt) = val;
        end
        function obj = set.d2(obj,val)
            obj.d2stor(:,1:obj.npt) = val;
        end
        function obj = set.h(obj,val)
            obj.h = val;
        end
        function dim = get.dim(obj)
            dim = size(obj.r,1);
        end
        function nptstor = get.nptstor(obj)
            nptstor = size(obj.rstor,2);
        end
        function datadim = get.datadim(obj)
            datadim = size(obj.data,1);
        end
        function npt = get.npt(obj)
            npt = obj.npt;
        end
        
        function obj = addpt(obj,nptadd)
            if nargin < 2
                nptadd = 1;
            end
            assert(and(isnumeric(nptadd),nptadd > 0), ...
                'nptadd must be positive integer');
            assert(obj.npt+nptadd <= obj.nptmax, ...
                'adding points would exceed maximum storage length');
            while obj.npt + nptadd > obj.nptstor
                % double size until sufficient storage available
                obj = resize(obj,max(min(2*obj.nptstor,obj.nptmax),...
                    min(obj.npt+nptadd,obj.nptmax)));
            end
            obj.npt = obj.npt+nptadd;
            assert((obj.npt <= min(obj.nptstor,obj.nptmax)),'something went wrong');
        end
            
        function obj = resize(obj,nptstornew)
            assert(nptstornew >= obj.npt, ...
                ['error new length is less than number of points\n', ...
                'info would be lost.']);
            assert(nptstornew <= obj.nptmax, ...
                ['new storage exceeds maximum storage.\n', ...
                'perhaps something went wrong. otherwise, increase nptmax']);
            rtemp = obj.r;
            dtemp = obj.d;
            d2temp = obj.d2;
            datatemp = obj.data;
            obj.rstor = zeros(obj.dim,nptstornew);
            obj.dstor = zeros(obj.dim,nptstornew);
            obj.d2stor = zeros(obj.dim,nptstornew);
            obj.datastor = zeros(obj.datadim,nptstornew);
            obj.r = rtemp;
            obj.d = dtemp;
            obj.d2 = d2temp;
            obj.data = datatemp;            
            obj.nptstor = nptstornew;
        end
        
        function obj = makedatarows(obj,nrows)
            if (nrows > 0)
                datatemp = obj.data;
                datadimold = obj.datadim;
                obj.datastor = zeros(datadimold+nrows,obj.nptstor);
                obj.data(1:datadimold,:) = datatemp(:,:);
                obj.hasdata = true;
            else
                if (nrows < 0)
                    warning('attempted to add negative rows, doing nothing');
                end
            end
        end
        
        function obj = cleardata(obj)
            obj.hasdata = false;
            obj.datastor = [];
        end
        
        [rn,dn,d2n,dist,tn,ichn] = nearest(obj,ref,ich,x,u)
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
        a = area(obj)
        s = arclength(obj)
    end
    methods(Static)
        obj = chunkfunc(fcurve,varargin)
    end
end
