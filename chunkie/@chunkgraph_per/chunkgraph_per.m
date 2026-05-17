classdef chunkgraph_per < chunkgraph

    %added properties: vstruc_use, per_dir
   properties
        vstruc_free
        merge_idx
   end
    methods

        function obj = chunkgraph_per(verts, edgesendverts, merge_idx, varargin)
            obj = obj@chunkgraph(verts, edgesendverts, varargin(:)); 
            obj = build_vstruc(obj,merge_idx); 
        end

        %note: multiple vert merging not supported
        function obj = build_vstruc(obj,merge_idx)
            obj.merge_idx = merge_idx; 
            obj.vstruc_free = obj.vstruc; 
            vstruc = obj.vstruc; 
            N_base_v = numel(vstruc); 
            N_merge = numel(merge_idx); 
            N_v_use = N_base_v - N_merge; 
            idx_skip = 1:N_base_v; 
             for i_merge = 1:N_merge
                        Nv_merge = numel(merge_idx{i_merge});
                        vert_merge = []; 
                        conn_merge = []; 
                        vert = merge_idx{i_merge}(1); 
                        edges = obj.edgesendverts; 
                        idx = edges(:) == merge_idx{i_merge};  
                        for i_vert = 1:Nv_merge
                            obj.edgesendverts(idx(:,i_vert)) = vert;
                            vert_merge = [vert_merge, vstruc{merge_idx{i_merge}(i_vert)}{1}]; 
                            conn_merge = [conn_merge, vstruc{merge_idx{i_merge}(i_vert)}{2}]; 
                        end
                        v_use{i_merge} = [{vert_merge} {conn_merge}]; 
                        idx_skip = setdiff(idx_skip,merge_idx{i_merge}); 
             end
             vrem = vstruc(idx_skip); 
             v_use = [v_use, vrem]; 
             obj.vstruc = v_use; 
    
             obj.v2emat = build_v2emat(obj); 
             obj.regions = findregions(obj); 
             obj = balance(obj); 
        end
    end
end