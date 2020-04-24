function [inds,adjs,info] = sortinfo(chnkr)
%SORTINFO computes a permutation of the chunk order in chnkr so that 
% adjacent chunks are consecutive. checks adjacency vector for errors.
% Provides diagnostics re: number of distinct curves in chunker,
% number of chunks per component, whether or not components are closed,
% and if there is an error in the adjacency info. Does not actually perform
% sort of the chunker
%
% Syntax: [inds,adjs,info] = sortinfo(chnkr)
%
% Input: 
%   chnkr - chunker object
%
% Output:
%   inds - indices which sort the chunks in chnkr so that adjacent chunks
%           are consecutive
%   adjs - sorted adjacency info (handles left and right ends, etc)
%   info - adjacency info structure
%       info.ncomp - number of curve components detected
%       info.nchs - number of chunks on each component
%       info.ifclosed - whether/or not each detected component is closed
%       info.ier - error flag
%           ier = 1, bad adj info, different number of left and right ends
%           ier = 2, bad adj info, missed/doubled chunks found
%
% Examples:
%   [inds,info] = sortinfo(chnkr)
%   [~,info] = sortinfo(chnkr)
%   inds = sortinfo(chnkr)
%
% see also SORT

% author: Travis Askham (askhamwhat@gmail.com)

nch = chnkr.nch;
ncompmax = 1000;

ier = 0;

adj2 = chnkr.adj(2,:);
adj1 = chnkr.adj(1,:);

ilefts = find(adj1 < 1);
irights = find(adj2 < 1);

if length(ilefts) ~= length(irights)
    % number of left and right ends is incompatible
    ier = 1;
    info = [];
    info.ier = ier;
    info.ncomp = 0;
    info.nchs = [];
    info.ifclosed = [];
    inds = 1:nch;
    adjs = chnkr.adj;
    return;
end

ifclosed=true(ncompmax,1);

% initialize

inds = -ones(nch,1);
ihit = false(nch,1);
nchs = zeros(ncompmax,1);
comp = 0;

adjs = [nch, 1:(nch-1);2:nch, 1];

ifclosed(1:length(ilefts)) = false;


% do components starting at vertices first
% attempt being contiguous (though that's impossible for some geometries)

ileftneg = find(adj1 < 0);
inext = 1;

for i = 1:length(ileftneg)
    comp = comp+1;
    icurrent = ileftneg(i);
    adjs(1,inext) = adj1(1,icurrent); %should always be negative here
    for ii = inext:nch
        inds(ii) = icurrent;
        ihit(icurrent) = true;
        itemp = adj2(icurrent);
        if (itemp == 0)
            nchs(comp) = ii+1-inext;
            inext = ii+1;
            adjs(2,ii) = 0;
            break;
        elseif (itemp < 0)
            % vertex
            nchs(comp) = ii+1-inext;
            inext = ii+1;
            adjs(2,ii) = itemp;
            % see if there's a component coming out of this 
            % vertex. if it exists, swap it to front
            iii = find(adj1(ileftneg(i+1:end)) == itemp,1,'first');
            if ~isempty(iii)
                iiitemp = ileftneg(i+iii);
                ileftneg(i+iii) = ileftneg(i+1);
                ileftneg(i+1) = iiitemp;
            end
            break;
        else
            icurrent = itemp;
        end
    end
end

ileftzero = find(adj1 == 0);
% then do components starting from a free edge
for i = 1:length(ileftzero)
    comp = comp+1;
    icurrent = ileftzero(i);
    adjs(1,inext) = adj1(icurrent); %should always be zero here
    for ii = inext:nch
        inds(ii) = icurrent;
        ihit(icurrent) = true;
        itemp = adj2(icurrent);
        if (itemp == 0)
            nchs(comp) = ii+1-inext;
            inext = ii+1;
            adjs(2,ii) = 0;
            break;
        elseif (itemp < 0)
            % vertex
            nchs(comp) = ii+1-inext;
            inext = ii+1;
            adjs(2,ii) = itemp;
        else
            icurrent = itemp;
        end
    end
end

% do remaining, closed components
for i = 1:ncompmax
    if inext > nch
        break;
    end
    comp = comp+1;
    icurrent = find(~ihit,1,'first');
    istart = icurrent;
    for ii = inext:nch
        inds(ii) = icurrent;
        ihit(icurrent) = true;
        itemp = adj2(icurrent);   
        if (itemp == istart)
            adjs(2,ii) = inext;
            adjs(1,inext) = ii;
            nchs(comp) = ii+1-inext;
            inext = ii+1;
            break;
        else
            icurrent = itemp;
        end
    end
end

ifclosed = ifclosed(1:comp);
nchs = nchs(1:comp);

if any(inds < 1) || length(unique(inds)) ~= length(inds)
    % sorting inds are bad
    ier = 2;
end

info.ifclosed = ifclosed;
info.nchs = nchs;
info.ncomp = comp;
info.ier = ier;
