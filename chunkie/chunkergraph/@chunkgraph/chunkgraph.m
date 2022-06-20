classdef chunkgraph

    properties(SetAccess=public)
        verts
        edge2verts
        echnks
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
    end   
    methods(Static)
        
    end
end
