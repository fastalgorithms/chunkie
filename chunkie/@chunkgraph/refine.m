function cg = refine(cg,opts,last_len)
%REFINE refine each edge of a chunkgraph object. 
%
% Syntax 
%   cg = refine(cg,opts);
%
% Input:
%   cg - chunkgraph
%   opts - cell array of structs, options to pass to chunker refine routine
%       in each face
%   last_len - optional float, ensure that the arclength of all panels 
%     touching each vertex is 2^-n * last_len for some n >= 0.
%
% Output:
%   cg - refined chunker, edges are sorted and balanced again
%
% see also CHUNKER/REFINE

if nargin == 1
    opts = [];
end

if length(opts) <= 1
    opts0 = opts;
    opts = cell(length(cg.echnks),1);
    for j = 1:length(cg.echnks)
        opts{j} = opts0;
    end
end

for j = 1:length(cg.echnks)
    chnkr = cg.echnks(j);
    chnkr = chnkr.refine(opts{j});
    chnkr = chnkr.sort();
    cg.echnks(j) = chnkr;
end

cg = balance(cg);


if nargin == 3
    % if a last length is provided
    nedges = length(cg.echnks);
    nverts = length(cg.vstruc);

    % first reparameterize each edge by arclength
    for i = 1:nedges
        len = sum(cg.echnks(i).wts(:));
        len_l = sum(cg.echnks(i).wts(:,1));
        len_r = sum(cg.echnks(i).wts(:,end));

        lvl_l = get_lvl(len_l, last_len);
        lvl_r = get_lvl(len_r, last_len);

        len_l = 2^-lvl_l*last_len;
        len_r = 2^-lvl_r*last_len;

        tsplit= [2*len_l, len - 2*len_r];
        tsplit = uniquetol(tsplit);

        mxlgth = max(sum(cg.echnks(i).wts(:,:)));
        cparams = []; cparams.maxchunklen = mxlgth;
        cparams.ta = 0; cparams.tb = len;
        cparams.tsplits = tsplit;
        [~,~,info] = sortinfo(cg.echnks(i)); 
        cparams.ifclosed = info.ifclosed;
        cparams.ifrefine = 1; cparams.lvlr = 't';

        param_data = chunkerarcparam_init(cg.echnks(i));
        fcurve = @(s) chunkerarcparam(s,param_data);

        cparams.eps = max(20 * param_data.eps, 1e-12);
        cg.echnks(i) = chunkerfunc(fcurve,cparams);
    end

    % now loop over vertices ensuring agreement
    for j = 1:nverts
        loc_edges = cg.vstruc{j}{1};
        loc_dir = cg.vstruc{j}{2};

        nloc = length(loc_edges);
        if nloc == 0
            continue
        end
        arcs = zeros(1,nloc);

        idch_loc = [cg.echnks(loc_edges).nch];
        idch_loc(loc_dir == -1) = 1;

        for k = 1:nloc
            wts = cg.echnks(loc_edges(k)).wts;
            arcs(k) = sum(wts(:,idch_loc(k)));
        end

        % find desired last panel size
        lvl = get_lvl(arcs, last_len);
        len = last_len * 2^-lvl;

        % find splitting fractions
        fracs = len./arcs;
        fracs(loc_dir==1) = 1 - fracs(loc_dir==1);
        if any(fracs<1e-8 | fracs>1-1e-8)
            if all(fracs<1e-8 | fracs>1-1e-8)
                % nothing to do
                continue;
            else
                % won't be able to find splits at this lvl, so go smaller
                lvl = lvl + 1;
                len = last_len * 2^-lvl;
                fracs = len./arcs;
                fracs(loc_dir==1) = 1 - fracs(loc_dir==1);
            end
        end

        for k = 1:nloc
            split_opt = []; 
            split_opt.frac = fracs(k);
            if split_opt.frac < 1e-8 || split_opt.frac > 1-1e-8
                continue, 
            end

            % split last panel to a multiple of last_len
            cg.echnks(loc_edges(k)) = split(cg.echnks(loc_edges(k)),idch_loc(k),split_opt);
            cg.echnks(loc_edges(k)) = sort(cg.echnks(loc_edges(k)));

            % now ensure that the last two panels are of correct length
            if loc_dir(k) == 1
                cg.echnks(loc_edges(k)) = split(cg.echnks(loc_edges(k)),cg.echnks(loc_edges(k)).nch);
            else
                cg.echnks(loc_edges(k)) = split(cg.echnks(loc_edges(k)),1);
            end
            cg.echnks(loc_edges(k)) = sort(cg.echnks(loc_edges(k)));
        end
    end

    % now that vertex panels are correct, we fix the rest
    cg = refine(cg);
end

end

function lvl = get_lvl(lens, last_len)
    % helper routine for finding the multiple of last_len to use in length
    % matcher
    lvl = -log2(min(lens)/last_len);
    if (lvl-round(lvl)) < 1e-8
        lvl = round(lvl);
    else
        lvl = ceil(lvl);
    end
    lvl = max(lvl, 1);
end
