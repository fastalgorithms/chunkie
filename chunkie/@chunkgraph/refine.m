function cg = refine(cg,opts)
%REFINE refine each edge of a chunkgraph object. 
%
% Syntax 
%   cg = refine(cg,opts);
%
% Input:
%   cg - chunkgraph
%   opts - options structure (default values)
%       opts.dlist = list of edges to refine (1:nedges)
%       opts.ilist = list of edges not to refine ([]), this overrules
%               opts.dlist
%       opts.last_len = if provided, ensure that the arclength of all 
%               panels touching each vertex is 2^-n * last_len for some 
%               n >= 0.
%       opts.splitchunks = cell array of integer arrays ([]), each element 
%               is a list of chunks to split, helpful for refining chunks 
%               where a function is not resolved
%       opts.nchmax = maximum number of chunks on refined chunker
%                   (10*chnkr.nch)
%       opts.lvlr = level restriction flag ('a'), if 'a', enforce that no
%                   two adjacent chunks should differ in length by more 
%                   than a factor of approx 2.1 (see lvlrfac). if 'n'
%                   don't enforce (not recommended)
%       opts.lvlrfac = (2.1, see chunker.lvlrfacdefault) factor for
%                      enforcing level restriction
%       opts.maxchunklen = maximum chunk length (Inf). enforce that no
%                   chunk is larger than this maximum length
%       opts.nover = oversample boundary nover times (0)
%       opts.maxiter_lvlr = number of iterations allowed when attempting
%                           level restriction (1000)
%
% Output:
%   cg - refined chunker, edges are sorted and balanced again
%
% see also CHUNKER/REFINE

if nargin == 1
    opts = [];
end

if ~isfield(opts,'ilist')
    opts.ilist = [];
end

if ~isfield(opts,'dlist')
    opts.dlist = 1:length(cg.echnks);
end
opts.dlist = setdiff(opts.dlist,opts.ilist);

if ~isfield(opts,'splitchunks')
    opts.splitchunks = [];
end
if ~iscell(opts.splitchunks)
    splitchunks = opts.splitchunks;
    opts.splitchunks = cell(1,length(cg.echnks));
    for i = 1:length(cg.echnks)
        opts.splitchunks{i} = splitchunks;
    end
end

for j = opts.dlist
    optsj = opts;
    optsj.splitchunks = opts.splitchunks{j};
    chnkr = cg.echnks(j);
    chnkr = chnkr.refine(optsj);
    chnkr = chnkr.sort();
    cg.echnks(j) = chnkr;
end

cg = balance(cg);


if isfield(opts,'last_len')
    % if a last length is provided
    if isempty(opts.last_len), return, end

    last_len = opts.last_len;
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
