function chnkr = refine(chnkr,opts)
%REFINE refine the chunker object according to any of a number of
% rules. This routine takes the chunker object to be *the domain* 
% i.e. starting with an under-resolved chunker representation of some 
% domain, you will not obtain a better approximation of that domain
% by refining the chunker in this way.
%
% Syntax: chnkr = refine(chnkr,opts)
%
% Input:
%   chnkr - chunker object
% Optional input:
%   opts - options structure (default values)
%       opts.splitchunks = list of chunks to split ([]), helpful for 
%                           refining chunks where a function is not 
%                           resolved
%       opts.nchmax = maximum number of chunks on refined chunker
%                   (10*chnkr.nch)
%       opts.lvlr = level restriction flag ('a'), if 'a', enforce that no
%                   two adjacent chunks should differ in length by more 
%                   than a factor of approx 2 (see lvlrfac). if 't', 
%                   enforce that the 
%                   length in parameter space of two adjacent chunks
%                   differs by no more than a factor of approx 2. if 'n'
%                   don't enforce (not recommended)
%       opts.lvlrfac = (2.0) factor for enforcing level restriction
%       opts.maxchunklen = maximum chunk length (Inf). enforce that no
%                   chunk is larger than this maximum length
%       opts.nover = oversample boundary nover times (0)
%       opts.maxiter_lvlr = number of iterations allowed when attempting
%                           level restriction (1000)
%
% Output:
%   chnkr - modified chunker object
%
% Examples:
%   opts = []; opts.lvlr = 'a'; opts.nover = 1;
%   %this will enforce level restriction on the curve and oversample once
%   chnkr = refine(chnkr,opts);
%
% see also CHUNKERFUNC, CHUNKPOLY, SPLIT

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 2
    opts = [];
end

lvlr = 'a';
lvlrfac = 2.0;
maxchunklen = Inf;
nover = 0;
splitchunks = [];

maxiter_lvlr=1000;
maxiter_maxlen=1000;

nchmax = chnkr.nchmax;

if isfield(opts,'lvlr'); lvlr = opts.lvlr; end
if isfield(opts,'lvlrfac'); lvlrfac = opts.lvlrfac; end
if isfield(opts,'maxchunklen'); maxchunklen = opts.maxchunklen; end
if isfield(opts,'nchmax'); nchmax = opts.nchmax; end
if isfield(opts,'nover'); nover= opts.nover; end
if isfield(opts,'splitchunks'); splitchunks = opts.splitchunks; end


nch = chnkr.nch;
k = chnkr.k;
dim = chnkr.dim;
vert = chnkr.vert;

% compute lengths of chunks at start and update along the way

nchlen = min(2*nch,nchmax);
chunklens = zeros(nchlen,1);
ws = weights(chnkr);
chunklens(1:nch) = sum(ws,1);

% 

[x,w,u] = lege.exps(k);


% splitting type

stype = 'a';
if strcmpi(lvlr,'t')
    stype = 't';
end


% chunks told to split

for i = 1:length(splitchunks)
    ii = splitchunks(i);

% split chunk ii now, and recalculate nodes, d, etc
    if (chnkr.nch + 1 > nchmax)
        error('too many chunks')
    end

    chnkr = split(chnkr,ii,[],x,w,u,stype);

    % update chunklens 

    nch = chnkr.nch;

    if (nch > length(chunklens))
        chunklens = resizechunklens(chunklens,nchmax);
    end

    di = chnkr.dstor(:,:,ii);
    hi = chnkr.hstor(ii);
    dnch = chnkr.dstor(:,:,nch);
    hnch = chnkr.hstor(nch);

    dsdti = sqrt(sum(di.^2,1))*hi;
    dsdtnch = sqrt(sum(dnch.^2,1))*hnch;

    chunklens(ii) = dot(dsdti,w);
    chunklens(nch) = dot(dsdtnch,w);

end    


% maximum chunklength

if maxchunklen < Inf
    for ijk = 1:maxiter_maxlen

        nchold=chnkr.nch;
        ifdone=1;

        for i = 1:nchold

            rlself = chunklens(i);

    %       only check if self is sufficiently small

            if rlself > maxchunklen

    %       split chunk i now, and recalculate nodes, d, etc

                if (chnkr.nch + 1 > nchmax)
                    error('too many chunks')
                end

                chnkr = split(chnkr,i,[],x,w,u,stype);

                % update chunklens 
                
                nch = chnkr.nch;
                
                if (nch > length(chunklens))
                    chunklens = resizechunklens(chunklens,nchmax);
                end
                
                di = chnkr.dstor(:,:,i);
                hi = chnkr.hstor(i);
                dnch = chnkr.dstor(:,:,nch);
                hnch = chnkr.hstor(nch);

                dsdti = sqrt(sum(di.^2,1))*hi;
                dsdtnch = sqrt(sum(dnch.^2,1))*hnch;

                chunklens(i) = dot(dsdti,w);
                chunklens(nch) = dot(dsdtnch,w);

                ifdone=0;

            end
        end

        if (ifdone == 1)
            break;
        end

    end
end

% level restriction

if (strcmpi(lvlr,'a') || strcmpi(lvlr,'t'))
    for ijk = 1:maxiter_lvlr

        nchold=chnkr.nch;
        ifdone=1;

        for i = 1:nchold
            i1=chnkr.adj(1,i);
            i2=chnkr.adj(2,i);

            if (strcmpi(lvlr,'t'))
                rlself = chnkr.h(i);
                rl1 = rlself;
                rl2 = rlself;
                if (i1 > 0)
                    rl1 = chnkr.h(i1);
                end
                if (i2 > 0)
                    rl2 = chnkr.h(i2);
                end
                if (i1 < 0)
                    % meets at vertex (parameter space ref not recommended)
                    rl1 = min(chnkr.h(vert{-i1}));
                end
                if (i2 < 0)
                    rl2 = min(chnkr.h(vert{-i2}));
                end
                    
            else

                rlself = chunklens(i);

                rl1=rlself;
                rl2=rlself;

                if (i1 > 0)
                    rl1 = chunklens(i1);
                end
                if (i2 > 0)
                    rl2 = chunklens(i2);
                end
                if (i1 < 0)
                    rl1 = min(chunklens(vert{-i1}));
                end
                if (i2 < 0)
                    rl2 = min(chunklens(vert{-i2}));
                end
            end

    %       only check if self is larger than either of adjacent blocks,
    %       iterating a couple times will catch everything

            
            sc = lvlrfac;
            if (rlself > sc*rl1 || rlself > sc*rl2)

    %       split chunk i now, and recalculate nodes, d, etc

                if (chnkr.nch + 1 > nchmax)
                    error('too many chunks')
                end

                chnkr = split(chnkr,i,[],x,w,u,stype);

                % update chunklens 

                nch = chnkr.nch;
                
                if (nch > length(chunklens))
                    chunklens = resizechunklens(chunklens,nchmax);
                end
                
                di = chnkr.dstor(:,:,i);
                hi = chnkr.hstor(i);
                dnch = chnkr.dstor(:,:,nch);
                hnch = chnkr.hstor(nch);

                dsdti = sqrt(sum(di.^2,1))*hi;
                dsdtnch = sqrt(sum(dnch.^2,1))*hnch;

                chunklens(i) = dot(dsdti,w);
                chunklens(nch) = dot(dsdtnch,w);

                ifdone=0;

            end
        end

        if (ifdone == 1)
            break;
        end

    end
end


% oversample 

for ijk = 1:nover

    nchold=chnkr.nch;

    for i = 1:nchold

%       split chunk i now, and recalculate nodes, d, etc
        if (chnkr.nch + 1 > nchmax)
            error('too many chunks')
        end

        chnkr = split(chnkr,i,[],x,w,u,stype);

        % update chunklens 

        nch = chnkr.nch;
        
        if (nch > length(chunklens))
            chunklens = resizechunklens(chunklens,nchmax);
        end
                
        di = chnkr.dstor(:,:,i);
        hi = chnkr.hstor(i);
        dnch = chnkr.dstor(:,:,nch);
        hnch = chnkr.hstor(nch);

        dsdti = sqrt(sum(di.^2,1))*hi;
        dsdtnch = sqrt(sum(dnch.^2,1))*hnch;

        chunklens(i) = dot(dsdti,w);
        chunklens(nch) = dot(dsdtnch,w);

    end

end

chnkr.n = normals(chnkr);

end


function chunklens = resizechunklens(chunklens,nchmax)

nchlen = length(chunklens);
nchnew = min(2*nchlen,nchmax);
chtemp = chunklens;
chunklens = zeros(nchnew,1);
chunklens(1:nchlen) = chtemp;

end