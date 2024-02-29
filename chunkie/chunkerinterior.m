function [in] = chunkerinterior(chnkobj,ptsobj,opts)
%CHUNKERINTERIOR returns an array indicating whether each point specified 
% % by pts is inside the domain. Assumes the domain is closed.
%
% Syntax: in = chunkerinterior(chnkobj,pts,opts)
%         in = chunkerinterior(chnkobj,{x,y},opts) % meshgrid version
%
% Input:
%   chnkobj - chunker object or chunkgraph object describing geometry
%   ptsobj - object describing the target points, can be specified as
%       * (chnkr.dim,:) array of points to test
%       * {x,y} - length 2 cell array. the points checked then have the
%           coordinates of a mesh grid [xx,yy] = meshgrid(x,y)
%       * chunker object, in which case it uses chunker.r(:,:) 
%       * chunkgraph object, in which case it uses chunkgraph.r(:,:) 
%
% Optional input:
%   opts - options structure with entries:
%       opts.fmm = boolean, use FMM 
%       opts.flam = boolean, use FLAM routines
%  Note on the default behavior: 
%    by default it tries to use the fmm if it exists, if it doesn't
%    then unless explicitly set to false, it tries to use flam
%
% Output:
%   in - logical array, if in(i) is true, then pts(:,i) is inside the
%       domain or for a mesh grid [xx(i); yy(i)] is inside the domain.
%
% Examples:
%   chnkr = chunkerfunc(@(t) starfish(t));
%   pts = 2*randn(2,100);
%   in = chunkerinterior(chnkr,pts);
%

% author: Travis Askham (askhamwhat@gmail.com)

grid = false;


if nargin < 3
    opts = [];
end

% Assign appropriate object to chnkr
if class(chnkobj) == "chunker"
   chnkr = chnkobj;
elseif class(chnkobj) == "chunkgraph"
   chnkr = merge(chnkobj.echnks);
else
    msg = "Unsupported object in chunkerinterior";
    error(msg)
end

assert(chnkr.dim == 2,'interior only well-defined for 2D');

% Assign appropriate object to pts
if isa(ptsobj, "cell")
    assert(length(ptsobj)==2,'second input should be either 2xnpts array or length 2 cell array');
    x = ptsobj{1};
    y = ptsobj{2};
    grid = true;
    [xx,yy] = meshgrid(x,y);
    pts = [xx(:).'; yy(:).'];
elseif isa(ptsobj, "chunker")
    pts = ptsobj.r(:,:);
elseif isa(ptsobj, "chunkgraph")
    pts = ptsobj.r(:,:);
else
    pts = ptsobj;
end


usefmm = true;
if isfield(opts,'fmm')
    usefmm = opts.fmm;
end

useflam = true;
if isfield(opts,'flam')
    useflam = opts.flam;
end

usefmm_final = false;
useflam_final = false;

if usefmm
    s = which('fmm2d');
    if(~isempty(s))
        usefmm_final = true;
    else
        useflam_final = useflam;
    end
else
   useflam_final = useflam;
end

% use bernstein ellipses and rectangles to flag problematic points
%

eps_local = 1e-3;

rho = 1.2;
optsflag = [];  optsflag.rho = rho; optsflag.occ = 5;
if grid
    flag = flagnear_rectangle_grid(chnkr,x,y,optsflag);
else
    flag = flagnear_rectangle(chnkr,pts,optsflag);
end

npoly = chnkr.k;
nlegnew = chnk.ellipse_oversample(rho,npoly,eps_local);
nlegnew = max(nlegnew,chnkr.k);

[chnkr2] = upsample(chnkr,nlegnew);

icont = false;
if usefmm_final
   try
       wchnkr = chnkr2.wts;
       dens1_fmm = ones(chnkr2.k*chnkr2.nch,1).*wchnkr(:);
       pgt = 1;
       vals1 = chnk.lap2d.fmm(eps_local,chnkr2,pts,'d',dens1_fmm,pgt);
   catch
       fprintf('using fmm failed due to incompatible mex, try regenrating mex\n');
       useflam_final = useflam;
       icont = true;
   end
end

if ~usefmm_final || icont
    kernd = kernel('lap','d');
    dens1 = ones(chnkr2.k,chnkr2.nch);
    

    opdims = [1 1];

    if useflam_final
        xflam1 = chnkr.r(:,:);
        matfun = @(i,j) chnk.flam.kernbyindexr(i,j,pts,chnkr2,kernd,opdims);
        [pr,ptau,pw,pin] = chnk.flam.proxy_square_pts();

        pxyfun = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
            ctr,chnkr2,kernd,opdims,pr,ptau,pw,pin);
        F = ifmm(matfun,pts,xflam1,200,1e-6,pxyfun);
        vals1 = ifmm_mv(F,dens1(:),matfun);
    else
        optskerneval = []; optskerneval.usesmooth = 1;
        vals1 = chunkerkerneval(chnkr2,kernd,dens1,pts,optskerneval);
    end
end
in = abs(vals1+1) < abs(vals1);

% for points where the integral might be inaccurate:
% find close boundary point and check normal direction


nnzpt = sum(flag~=0,2);
ipt = find(nnzpt);

npts = numel(pts)/2;
distmins = inf(npts,1);
dss = zeros(2,npts);
rss = zeros(2,npts);

k = chnkr.k;
[t,~,u] = lege.exps(k);

for i = 1:chnkr.nch

    % check side based on closest boundary node
    rval = chnkr.r(:,:,i);
    dval = chnkr.d(:,:,i);
    nval = chnkr.n(:,:,i);
    [ji] = find(flag(:,i));
    ptsi = pts(:,ji);
    nptsi = size(ptsi,2);
    dist2all = reshape(sum( abs(reshape(ptsi,2,1,nptsi) ...
                    - reshape(rval,2,k,1)).^2, 1),k,nptsi);
    [dist2all,ipti] = min(dist2all,[],1);
    for j = 1:length(ji)
        jj = ji(j);
        if dist2all(j) < distmins(jj)
            distmins(jj) = dist2all(j);
            rss(:,jj) = rval(:,ipti(j));
            dss(:,jj) = dval(:,ipti(j));
        end
    end

    % if angle is small do a refined search for closest point
    ptsn = nval(:,ipti);
    diffs = rval(:,ipti)-ptsi;
    dots = sum(ptsn.*diffs,1);
    
    jsus = abs(dots) < 2e-1*sqrt(dist2all);
    jii = ji(jsus);
    [~,rs,ds,~,dist2s] = chnk.chunk_nearparam(rval,pts(:,jii),[],t,u);
    for j = 1:length(jii)
        jj = jii(j);
        if dist2s(j) < distmins(jj)
            distmins(jj) = dist2s(j);
            rss(:,jj) = rs(:,j);
            dss(:,jj) = ds(:,j);
        end
    end
end

for i = 1:length(ipt)
    jj = ipt(i);
    in(jj) = (rss(:,jj)-pts(:,jj)).'*[dss(2,jj);-dss(1,jj)] > 0;
end
