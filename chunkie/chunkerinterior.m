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
%       opts.axissym = boolean, chunker is axissymmetric
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
elseif isa(ptsobj, "chunker") || isa(ptsobj, "chunkgraph") || ...
        (isstruct(ptsobj) && isfield(ptsobj,"r"))
    pts = ptsobj.r(:,:);
elseif isnumeric(ptsobj)
    pts = ptsobj;
else
    msg = "chunkerinterior: input 2 not a recognized type";
    error(msg);
end


usefmm = true;
if isfield(opts,'fmm')
    usefmm = opts.fmm;
end

useflam = true;
if isfield(opts,'flam')
    useflam = opts.flam;
end

axissym = false;
if isfield(opts,'axissym')
    axissym = opts.axissym;
end

if axissym
    nch = chnkr.nch;
    istart = nch+1;
    iend = 2*nch;
    chnkr = sort(chnkr);
    chnkr = chnkr.addchunk(nch);
    chnkr.r(:,:,istart:iend) = fliplr(chnkr.r(:,:,1:nch));
    chnkr.r(1,:,istart:iend) = -chnkr.r(1,:,istart:iend);
    chnkr.d(:,:,istart:iend) = fliplr(chnkr.d(:,:,1:nch));
    chnkr.d(2,:,istart:iend) = -chnkr.d(2,:,istart:iend);

    chnkr.d2(:,:,istart:iend) = fliplr(chnkr.d2(:,:,1:nch));
    chnkr.d2(1,:,istart:iend) = chnkr.d2(1,:,istart:iend);
    
    chnkr.wts = weights(chnkr);
    chnkr.n = normals(chnkr);
    chnkr.adj(1,:) = 0:chnkr.nch-1;
    chnkr.adj(2,:) = 2:chnkr.nch+1;
    chnkr.adj(1,1) = chnkr.nch;
    chnkr.adj(2,chnkr.nch) = 1;
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


eps_local = 1e-3;
rho = 1.6;
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

iffy = min(abs(vals1+1),abs(vals1)) > 1e-2;

ipt = find(iffy(:));
pts_iffy = pts(:,iffy);

optsflag = [];  optsflag.rho = rho; optsflag.occ = 5;
flag = flagnear_rectangle(chnkr,pts_iffy,opts);
flag2 = flagnear_rectangle(chnkr,chnkr.r(:,:),opts);
flag2 = ((flag2.'*kron(speye(chnkr.nch),ones(chnkr.k,1))) > 1).';

flag = (flag2*flag.').';

npts_iffy = numel(pts_iffy)/2;
assert(npts_iffy == length(ipt));
distmins = inf(npts_iffy,1);
dss = zeros(2,npts_iffy);
rss = zeros(2,npts_iffy);

k = chnkr.k;
[t,~,u] = lege.exps(k);

for i = 1:chnkr.nch

    % check side based on closest boundary node
    rval = chnkr.r(:,:,i);
    dval = chnkr.d(:,:,i);
    nval = chnkr.n(:,:,i);
    [ji] = find(flag(:,i));

    ptsi = pts_iffy(:,ji);
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
    
    jsus = abs(dots) < 0.9*sqrt(dist2all);
    jii = ji(jsus);
    [~,rs,ds,~,dist2s] = chnk.chunk_nearparam(rval,pts_iffy(:,jii),[],t,u);
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
    if distmins(i) < inf
        jj = ipt(i);
        in(jj) = (rss(:,i)-pts_iffy(:,i)).'*[dss(2,i);-dss(1,i)] > 0;
    end
end
