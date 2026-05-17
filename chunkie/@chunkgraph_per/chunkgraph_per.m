classdef chunkgraph_per < chunkgraph

    %added properties: vstruc_use, per_dir
   properties
    vstruc_use
   end
    methods

        function obj = chunkgraph_per(verts, edgesendverts, merge_idx, varargin)
            obj = obj@chunkgraph(verts, edgesendverts, varargin(:)); 
            obj = build_vstruc(obj,merge_idx); 
        end

        %note: multiple vert merging not supported
        function obj = build_vstruc(obj,merge_idx)
        vstruc = obj.vstruc; 
        N_base_v = numel(vstruc); 
        N_merge = numel(merge_idx); 
        N_v_use = N_base_v - N_merge; 
        v_use = cell(1,N_v_use); 
        idx_skip = 1:N_base_v; 
         for i_merge = 1:N_merge
                    Nv_merge = numel(merge_idx{i_merge}); 
                    vert_merge = zeros(1,Nv_merge); 
                    conn_merge = zeros(1,Nv_merge); 
                    for i_vert = 1:Nv_merge
                        vert_merge(i_vert) = vstruc{merge_idx{i_merge}(i_vert)}{1}; 
                        conn_merge(i_vert) = vstruc{merge_idx{i_merge}(i_vert)}{2}; 
                    end
                    v_use{i_merge} = [{vert_merge} {conn_merge}]; 
                    idx_skip = setdiff(idx_skip,merge_idx{i_merge}); 
         end
         v_use(N_merge+1:end) = vstruc(idx_skip);

         obj.vstruc_use = v_use; 
        end
    end
end