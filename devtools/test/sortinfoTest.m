sortinfoTest0();


function sortinfoTest0()
% test that sortinfo gives the same inds, adjs and info as the old
% sortinfo on a perforated disk (outer circle plus many small holes, so
% one closed component per hole)
%
% sortinfo_old, defined below, is the older version of sortinfo (before
% the fix), kept here so that the new sortinfo can be compared to it
%
% to time the old and new sortinfo, run this test with larger values of l
% (e.g. ls = 10:10:60) and put tic/toc around the two sortinfo calls

ls = 10:10:30;

for j = 1:length(ls)

    [chnkr, nholes] = get_chunker(ls(j));

    [inds1, adjs1, info1] = sortinfo_old(chnkr);
    [inds2, adjs2, info2] = sortinfo(chnkr);

    % checks the correctness of the code

    assert(isequal(inds1, inds2) && isequal(adjs1, adjs2) && isequal(info1, info2));
    assert(info2.ier == 0 && info2.ncomp == nholes + 1);

end

end


%% Old Sortinfo

function [inds,adjs,info] = sortinfo_old(chnkr)
% SORTINFO computes a permutation of the chunk order in chnkr so that 
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

ncompmax = nch;

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
    icurrent = find(~ihit, 1, 'first');
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
end


%% Definition of geometry


function [chnkr, nholes] = get_chunker(j)
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];



l = j;
k = 1/(4*l);
h = 5*l; 

rmax = 2/h;
nch = max(10, ceil(2*pi/rmax) + 2);

opts = [];
opts.maxchunklen = 1.1;
pref = []; pref.nchmax = 20000;

chnkr1 = chunkerfuncuni(circfun, nch, opts, pref);


chunkerhole0 = k*chunkerfuncuni(circfun, 6);   
chunkerhole0 = chunkerhole0.reverse();
nh = 0;
chunkers(10000,1) = chunker();
chunkers(1) = chnkr1;
ictr = 1;
for i = 0 : l-1
   for j = 0 : l-1
        check1 = [i/l,j/l];
        l1 = norm(check1 + [1/h,1/h]); 
        l2 = norm(check1 + [1/l,0.0] + [-1/h,1/h]);
        l3 = norm(check1 + [0.0,1/l] + [1/h,-1/h]);
        l4 = norm(check1 + [1/l,1/l] + [-1/h,-1/h]);
        if (l1<=1) && (l2<=1) && (l3<=1) && (l4<=1)

            chnkrctr1 = check1 + (1/(2*l))*[1.0,1.0];
            chnkrctr2 = [-1,1].*chnkrctr1;
            chnkrctr3 = [1,-1].*chnkrctr1;
            chnkrctr4 = [-1,-1].*chnkrctr1;
            chnkrhole1 = chnkrctr1 + chunkerhole0;
            chnkrhole2 = chnkrctr2 + chunkerhole0;
            chnkrhole3 = chnkrctr3 + chunkerhole0;
            chnkrhole4 = chnkrctr4 + chunkerhole0;
            chunkers(ictr+1) = chnkrhole1;
            chunkers(ictr+2) = chnkrhole2;
            chunkers(ictr+3) = chnkrhole3;
            chunkers(ictr+4) = chnkrhole4;
            ictr = ictr + 4;
            nh = nh + 4;
         end
    end
end
nholes = nh;
pref = [];
pref.nchmax = chnkr1.nch + 6*nh;

% merge the holes in blocks of 100

nblock = 100;
blocks = chunker();
iblock = 0;
for m0 = 1:nblock:nholes
    m1 = min(m0+nblock-1, nholes);
    holes = chunker();
    for m = m0:m1
        holes(m-m0+1,1) = chunkers(m+1);
    end
    iblock = iblock+1;
    blocks(iblock,1) = merge(holes, pref);
end

% outer circle first, then all the holes
chnkr = merge([chnkr1; blocks], pref);
end
