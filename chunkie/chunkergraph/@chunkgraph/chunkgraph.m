classdef chunkgraph

    properties(SetAccess=private)
        edge2verts
        verts
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
            obj.verts = val;
        end
    end   
    methods(Static)
        
    end
end
