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
    end   
    methods(Static)
        
    end
end
