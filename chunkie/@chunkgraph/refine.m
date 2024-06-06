function cg = refine(cg,opts)
%REFINE refine each edge of a chunkgraph object. 
%
% Syntax 
%   cg = refine(cg,opts);
%
% Input:
%   cg - chunkgraph
%   opts - struct, options to pass to chunker refine routine
%
% Output:
%   cg - refined chunker, edges are sorted and balanced again
%
% see also CHUNKER/REFINE

for j = 1:length(cg.echnks)
    chnkr = cg.echnks(j);
    chnkr = chnkr.refine(opts);
    chnkr = chnkr.sort();
    cg.echnks(j) = chnkr;
end

cg = balance(cg);