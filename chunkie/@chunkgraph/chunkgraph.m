classdef chunkgraph

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
        adj
        sourceinfo
        npts
    end
    
    methods
        function obj = chunkgraph(p)
        
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
        function adj = get.adj(obj)
            chnk = merge(obj.echnks);
            adj = chnk.adj;
        end
        function npts = get.npts(obj)
            npts = 0;
            for iedge = numel(obj.echnks)
            	n = size(obj.echnks(iedge).r(:,:),2);
                npts = n + npts;
            end
        end
        function sourceinfo = get.sourceinfo(obj)
            sourceinfo = [];
            ntot = obj.npts;
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
