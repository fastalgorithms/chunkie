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
% obj = calc_per(obj,merge_idx)
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
        vert_per %NOT CURRENTLY USED... probably will use in chunkermat
        dx
        dy
   end
    methods

        function obj = chunkgraph_per(verts, edgesendverts, merge_idx, varargin)
            obj = obj@chunkgraph(verts, edgesendverts, varargin{:}); 
            obj = calc_per(obj,merge_idx,varargin{:}); 
            obj.merge_idx = merge_idx; 
            if ~isempty(merge_idx)
                obj = build_vstruc(obj,merge_idx); 
            else
                obj.vstruc_free = obj.vstruc; 
                obj.edgesendverts_free = obj.edgesendverts; 
                obj.vert_per = nan(size(obj.verts)); 
            end
            obj.regions = findregions(obj); 
            obj = balance(obj); 
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
                vper = nan(size(obj.verts)); 
                nmerge = size(merge_idx,2); 
                dx = []; dy = []; 
                for m = 1:nmerge
                    vm = merge_idx{m}; 
                    v = obj.verts(:,vm); 
                    dxmat = abs(v(1,:)'-v(1,:)); 
                    dymat = abs(v(2,:)'-v(2,:));
                    dxm = unique(nonzeros(dxmat)); 
                    dym = unique(nonzeros(dymat)); 
                    if numel(dxm)>1 || numel(dym)>1
                        error('Periodicity of object is not consistent.')
                    end
                    if ~isempty(dxm)
                        vper(1,vm) = dxm; 
                    end
                    if ~isempty(dym)
                        vper(2,vm) = dym; 
                    end
                    dx = [dx;dxm]; 
                    dy = [dy;dym]; 
                end
                dx = unique(dx); dy = unique(dy); 
                if (numel(dx)>1) || (numel(dy)>1)
                    error('Periodicity of object is not consistent.')
                end
                obj.vert_per = vper; 
                obj.dx = dx; obj.dy = dy; 

            elseif isfield(cparams,'dx') && isfield(cparams,'dy')
                   obj.dx = cparams.dx; obj.dy = cparams.dy; 
            elseif isfield(cparams,'dx')
                obj.dx = cparams.dx; obj.dy = 0; 
            elseif isfield(cparams,'dy')
                obj.dx = 0; obj.dy = cparams.dy; 
            else
                error('Please provide merge_idx or period of closed geometry')
            end


        end

        function obj = build_vstruc(obj,merge_idx)
            obj.vstruc_free = obj.vstruc; 
            obj.edgesendverts_free = obj.edgesendverts;
            vstruc = obj.vstruc; 
            Nbase_v = numel(vstruc); 
            Nmerge = numel(merge_idx); 
            idx_skip = 1:Nbase_v; 
             for m = 1:Nmerge
                        Nv_merge = numel(merge_idx{m});
                        vert_merge = []; 
                        conn_merge = []; 
                        basevert = merge_idx{m}(1); 
                        idx = obj.edgesendverts(:) == merge_idx{m};  
                        for i_vert = 1:Nv_merge
                            obj.edgesendverts(idx(:,i_vert)) = basevert;
                            vert_merge = [vert_merge, vstruc{merge_idx{m}(i_vert)}{1}]; 
                            conn_merge = [conn_merge, vstruc{merge_idx{m}(i_vert)}{2}]; 
                        end
                        v_use{m} = [{vert_merge} {conn_merge}]; 
                        idx_skip = setdiff(idx_skip,merge_idx{m}); 
             end
             vrem = vstruc(idx_skip); 
             v_use = [v_use, vrem]; 
             obj.vstruc = v_use; 
             obj.v2emat = build_v2emat(obj); 
        end
    end
end