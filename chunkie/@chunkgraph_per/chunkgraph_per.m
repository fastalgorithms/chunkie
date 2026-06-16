classdef chunkgraph_per < chunkgraph
%CHUNKGRAPH_PER subclass of CHUNKGRAPH for geometries that are periodic in
% one or two directions. The period cell is described by the same vertex/
% edge graph as CHUNKGRAPH, but pairs of boundary vertices are identified
% to encode the periodicity. See CHUNKGRAPH for a full description
% of inherited properties and methods.
%
% The key addition is the merge_idx argument, which specifies a list
% vertices to be identified. The constructor computes the x and y periods
% automatically.
%
% chunkgraph_per properties (in addition to those of CHUNKGRAPH):
%   dx            - period in the x direction (0 if not periodic in x)
%   dy            - period in the y direction (0 if not periodic in y)
%   merge_idx     - (1 x nmerge) cell array; merge_idx{m} is a vector of
%                     vertex indices that are identified with each other
%   vstruc_free   - unmerged copy of vstruc
%   edgesendverts_free - unmerged copy ofedgesendverts
%   vert_per      - (2 x nverts) array; vert_per(:,j) is the displacement
%                     from vertex j to the first vertex of its merge group
%                     (merge_idx{m}(1)), or NaN if vertex j is not merged
%
% chunkgraph_per methods (in addition to those of CHUNKGRAPH):
%   obj = chunkgraph_per(verts,edgesendverts,merge_idx,varargin)
%   obj = findregions(obj) - region labeling aware of periodic identification
%   plot_regions(obj,iflabel) - plot with period-cell boundary shown
%   cg  = merge(cgs) - merge an array of chunkgraph_per objects into one
%   cg  = stack_layers(cg,merge_idx) - stack chunkgraph_per layers vertically
%   [dom,merge_idx] = gen_comp_domain(obj,varargin) - generate computational
%                       domain (scatter geometry + unit cell boundary)
%   [unbnd,bnd] = classify_loops(obj,loops) - split loops into unbounded
%                       periodic curves and bounded/closed objects
%   d   = loop_displacement(obj,edges) - net displacement around an edge loop
%   m   = loop_max_jump(obj,edges) - largest gap between consecutive edges
%   y   = loop_mean_y(obj,edges) - mean y-coordinate along an edge list
%   ny  = loop_normal_y(obj,edges) - mean outward normal y-component
%   poly = cell_polygon(obj,edges) - closed polygon for a tiling-closed loop
%
% Syntax:
%
%   cgp = chunkgraph_per(verts, edgesendverts, merge_idx)
%   cgp = chunkgraph_per(verts, edgesendverts, merge_idx, fchnks, cparams, pref)
%   cgp = chunkgraph_per(verts, edgesendverts, [], fchnks, cparams)
%
% Input:
%   verts          - (2 x nverts) vertex locations; same format as CHUNKGRAPH
%   edgesendverts  - (2 x nedges) starting/ending vertex indices per edge
%   merge_idx      - (1 x nmerge) cell array of vertex groups to identify,
%                      or [] if the period is instead specified via cparams
% Optional input:
%   fchnks, cparams, pref - same as CHUNKGRAPH
%   cparams.dx, cparams.dy - period(s), used when merge_idx is empty and
%                      the geometry consists of closed loops inside the cell
%
% see also CHUNKGRAPH, CHUNKGRAPH_PERINREGION, FINDREGIONS



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
                    dxmat = v(1,:)'-v(1,:); 
                    dymat = v(2,:)'-v(2,:);
                    dxm = unique(abs(nonzeros(dxmat))); 
                    dym = unique(abs(nonzeros(dymat))); 
                    if numel(dxm)>1 || numel(dym)>1
                        error('Periodicity of object is not consistent.')
                    end
                    if ~isempty(dxm)
                        vper(1,vm) = nonzeros(dxmat.'); 
                    end
                    if ~isempty(dym)
                        vper(2,vm) = nonzeros(dymat.'); 
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
            % merge chunkgraph vstrucs using periodicity
            obj.vstruc_free = obj.vstruc; 
            obj.edgesendverts_free = obj.edgesendverts;
            vstruc = obj.vstruc;
            for m = 1:numel(merge_idx)
                vm = merge_idx{m};
                basevert = vm(1);
            
                vert_merge = [];
                conn_merge = [];
            
                for ii = 1:numel(vm)
                    vcur = vm(ii);
            
                    vert_merge = [vert_merge, vstruc{vcur}{1}];
                    conn_merge = [conn_merge, vstruc{vcur}{2}];
            
                    idx = obj.edgesendverts == vcur;
                    obj.edgesendverts(idx) = basevert;
                end
            
                vstruc{basevert} = {vert_merge, conn_merge};
            
                for ii = 2:numel(vm)
                    vstruc{vm(ii)} = [];
                end
            end
            
            obj.vstruc = vstruc;
            obj.v2emat = build_v2emat(obj);
        end
    end
end