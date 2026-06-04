classdef chunkgraph_per < chunkgraph
% CHUNKGRAPH_PER is a subclass of CHUNKGRAPH that extends periodically in
% at least one direction. 
%
% INPUTS: same as chunkgraph and additionally: 
% -merge_idx: cell array with vecs, listing vertices to be merged
%
% properties: 
% -vstruc_free: original unmerged vstruc from chunkgraph
% -merge_idx: cell array, each cell containing vec of vertices to be merged
% -vert_per: vec housing periods associated with each vertex (each belonging to
% a periodic region). For vertices not being merged, vert_per(vert) = NaN.
%
% methods: 
% obj = chunkgraph_per(verts, edgesendverts, merge_idx, varargin)
%       -initializations chunkgraph_per object
% obj = build_vstruc(obj,merge_idx); 
%       - computes vert_per for each vertex
%       - loops over each cell in merge_idx, merges indices and overwrites
%       vstruc
%
%OUTPUT: 



   properties
        vstruc_free
        edgesendverts_free
        merge_idx
        vert_per
        dx
        dy
   end
    methods

        function obj = chunkgraph_per(verts, edgesendverts, merge_idx, varargin)

            obj = obj@chunkgraph(verts, edgesendverts, varargin{:}); 
            obj = calc_per(obj,merge_idx,varargin{:}); 
            if ~isempty(merge_idx)
                obj = build_vstruc(obj,merge_idx); 
            else
                obj.vstruc_free = obj.vstruc; 
                obj.edgesendverts_free = obj.edgesendverts; 
                obj.vert_per = nan(1,length(obj.verts(1,:))); 
            end
        end

        function obj = calc_per(obj,merge_idx,varargin)
         
            if ~isempty(varargin)
                if length(varargin)>1
                    cparams = varargin{2}; 
                end
            else
                cparams = []; 
            end
            if ~isempty(merge_idx)
                nmerge = size(merge_idx,2); 
                dx = []; dy = []; 
                for m = 1:nmerge
                    vm = merge_idx{m}; 
                    v = obj.verts(:,vm); 
                    dxmat = abs(v(1,:)'-v(1,:)); 
                    dymat = abs(v(2,:)'-v(2,:));
                    dx = [dx;unique(nonzeros(dxmat))]; 
                    dy = [dy;unique(nonzeros(dymat))]; 
                end
                dx = unique(dx); dy = unique(dy); 
                if (numel(dx)>1) || (numel(dy)>1)
                    error('Periodicity of object is not consistent.')
                end
                obj.dx = dx; obj.dy = dy; 

            elseif isfield(cparams,'dx') && isfield(cparams,'dy')
                   obj.dx = cparams.dx; obj.dy = cparams.dy; 
            elseif isfield(cparams,'dx')
                obj.dx = cparams.dx; obj.dy = []; 
            elseif isfield(cparams,'dy')
                obj.dx = []; obj.dy = cparams.dy; 
            else
                error('Please provide merge_idx or period of closed geometry')
            end


        end

        function obj = build_vstruc(obj,merge_idx)
            obj.vstruc_free = obj.vstruc; 
            obj.edgesendverts_free = obj.edgesendverts;
            vper = nan(numel(obj.verts(1,:)),1); 

            obj.merge_idx = merge_idx; 
            vstruc = obj.vstruc; 
            N_base_v = numel(vstruc); 
            N_merge = numel(merge_idx); 
            idx_skip = 1:N_base_v; 
             for i_merge = 1:N_merge
                        Nv_merge = numel(merge_idx{i_merge});

                        vert_merge = []; 
                        conn_merge = []; 
                        mi = merge_idx{i_merge}; 
                        verts = obj.verts(:,mi);
                        d = sqrt((verts(1,1)-verts(1,2))^2 + (verts(2,1)-verts(2,2))^2);
                        vper(mi) = d; 
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

             obj.vert_per = vper; 
             obj.vstruc = v_use;
            
             obj.v2emat = build_v2emat(obj); 
             obj0 = obj; obj0.vstruc = obj.vstruc_free; 
             obj0.edgesendverts = obj.edgesendverts_free; 
             obj.regions = findregions(obj); 
             obj = balance(obj); 
        end
    end
end