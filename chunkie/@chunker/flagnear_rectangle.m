function flag = flagnear_rectangle(chnkr,pts,opts)
%FLAGNEAR_RECTANGLE flag points which require special quadrature
% (or given factor of a chunklength) of any point on each chunk. on return is a sparse
% logical array. If the (i,j) entry is non zero then pts(:,i) is close
% to chunk j in chnkr.
%
% Syntax: flag = flagnear(chnkr,pts,opts)
%
% Input:
%   chnkr - chunker object describing curve
%   pts - (chnkr.dim,nt) array of point coordinates
%
% Optional input:
%   opts - options structure
%       opts.rho = Bernstein ellipse parameter (default=1.8)
%
% Output:
%   flag - (nt,chnkr.nch) sparse array. a non zero entry (i,j) means 
%       that the distance from pts(:,i) to at least one node on 
%       chunk j is less than opts.fac*length of chunk j.
%

% TODO: upgrade to hung chunk tree in current form this can be 
% slow for coarser boundaries
% author: Travis Askham (askhawhat@gmail.com)

% find nearest neighbors at certain level of refinement (here chosen
% uniformly, for simplicity)

rho = 1.8;
occ = 5;

if nargin < 3
    opts = [];
end

if isfield(opts,'rho')
    rho = opts.rho;
end
if isfield(opts,'occ')
    occ = opts.occ;
end

npt = chnkr.npt;
nch = chnkr.nch;
k = chnkr.k;
dim = chnkr.dim;

% the regions needing special quadrature are images of 
% bernstein ellipses

ells = ellipses(chnkr,rho);

[~,~,u] = lege.exps(k);
p0 = lege.pols(0,k-1);
val0 = (p0(:)).'*u;
dperm = permute(chnkr.d,[2,1,3]); dperm = dperm(:,:);
d0 = val0*dperm;
d0 = reshape(d0,chnkr.dim,chnkr.nch);
d0nrm = sqrt(sum(d0.^2,1));
d0 = d0./d0nrm;
d1s = d0;
d2s = flipud(d1s);
d2s(2,:) = -d2s(2,:);

%[rects,rectinfo] = bounding_rects(ells);
[rects,rectinfo] = bounding_rects_cheap(ells,d1s,d2s);

% make a "bounding volume hierarchy" of the regions
% needing special quadrature.
% T - a tree built on chunk centroids
% - For each leaf, build a grid aligned rectangle containing
%   the special quadrature regions
% - For nodes with children, find the rectangle containing the rectangles
%   of the box's children

ctrs = centroids(chnkr);
T = hypoct(ctrs,occ);

nnodes = length(T.nodes);
nlev = T.nlvl;
bvhbounds = zeros(dim,2,nnodes);
bvhbounds(:,2,:) = -1;


for l = nlev:-1:1
    inds = T.lvp(l)+1:T.lvp(l+1);
    for i = 1:length(inds)
        ii = inds(i);
        if isempty(T.nodes(ii).chld)
            if ~isempty(T.nodes(ii).xi)
                ellpts = rects(:,:,T.nodes(ii).xi);
                bvhbounds(:,1,ii) = min(ellpts(:,:),[],2);
                bvhbounds(:,2,ii) = max(ellpts(:,:),[],2);
            end
        else
            bvpts = bvhbounds(:,:,T.nodes(ii).chld);
            bvhbounds(:,1,ii) = min(bvpts(:,:),[],2);
            bvhbounds(:,2,ii) = max(bvpts(:,:),[],2);
        end
    end
end


nnzero = 0;
nn = 3*length(pts);
isp = zeros(nn,1);
jsp = zeros(nn,1);
istack = zeros(4*nlev,1);

for i = 1:size(pts,2)
    ntry = 0;
    is = 1;
    pt = pts(:,i);
    x = pt(1); y = pt(2);
    istack(1)=1;
    while(and(is > 0,ntry <= nnodes))
        ntry = ntry + 1;
        inode = istack(is);
        bvhtmp = bvhbounds(:,:,inode);
        xl = bvhtmp(1,1); xu = bvhtmp(1,2);
        yl = bvhtmp(2,1); yu = bvhtmp(2,2);
        if (x >= xl && x <= xu && y >= yl && y <= yu)
            %inside this box
            chld = T.nodes(inode).chld;
            if ~isempty(chld)
                % if has children, add them to stack
                istack(is:is+length(chld)-1) = chld;
                is = is + length(chld)-1;
            else
                xi = T.nodes(inode).xi;
                if ~isempty(xi)
                    it = [];
                    js = [];
                    for jj = 1:length(xi)
                        jell = xi(jj);
                        %in = inpolygon(pt(1),pt(2), ...
                        %     ells(1,:,jell),ells(2,:,jell));
                        d1 = pt(1)*rectinfo(1,1,jell)+ ...
                            pt(2)*rectinfo(2,1,jell);
                        d2 = pt(1)*rectinfo(1,2,jell)+ ...
                            pt(2)*rectinfo(2,2,jell);
                        in = (d1 >= rectinfo(1,3,jell) && ...
                            d1 <= rectinfo(2,3,jell) && ...
                            d2 >= rectinfo(1,4,jell) && ...
                            d2 <= rectinfo(2,4,jell));
                        if (in)
                            it = [it; i];
                            js = [js; jell];
                        end
                    end
                    nnew = length(it);
                    if nnew + nnzero > nn
                        itemp = isp;
                        jtemp = jsp;
                        isp = zeros(2*nn,1);
                        jsp = zeros(2*nn,1);
                        isp(1:nn) = itemp;
                        jsp(1:nn) = jtemp;
                        nn = 2*nn;
                    end
                    isp(nnzero+1:nnzero+nnew) = it;
                    jsp(nnzero+1:nnzero+nnew) = js;
                    nnzero = nnew + nnzero;
                    
                end
                is = is-1; %move up stack
            end
        else
            is = is-1; % move up stack
        end
    end
end

vsp = ones(nnzero,1);
isp = isp(1:nnzero);
jsp = jsp(1:nnzero);
flag = sparse(isp,jsp,vsp,length(pts),nch);

end

function [rects,rectinfo] = bounding_rects(convreg)
% find minimal area bounding rectangle for each region

[dim,m,n] = size(convreg);
xreg = reshape(convreg(1,:,:),m,n);
yreg = reshape(convreg(2,:,:),m,n);

dx = diff([xreg; xreg(1,:)],1,1);
dy = diff([yreg; yreg(1,:)],1,1);
dnrm = sqrt(dx.^2+dy.^2);

dx = reshape(dx./dnrm,1,m,n);
dy = reshape(dy./dnrm,1,m,n);

d1 = [dx;dy];
d2 = [-dy;dx];

rects = zeros(2,4,n);
rectinfo = zeros(2,4,n);

for i = 1:n
    d1i = d1(:,:,i);
    d2i = d2(:,:,i);
    pts = convreg(:,:,i).';
    d1c = pts*d1i;
    d2c = pts*d2i;
    d1emax = max(d1c,[],1);
    d1emin = min(d1c,[],1);
    d1eh = d1emax-d1emin;
    d2emax = max(d2c,[],1);
    d2emin = min(d2c,[],1);
    d2eh = d2emax-d2emin;
    areas = reshape(d1eh.*d2eh,m,1);
    [~,j] = min(areas(:));
    
    d1jmax = d1emax(j);
    d1jmin = d1emin(j);
    d2jmax = d2emax(j);
    d2jmin = d2emin(j);
    
    d1j = d1i(:,j);
    d2j = d2i(:,j);
    
    rects(:,:,i) = [d1jmax*d1j+d2jmax*d2j, d1jmin*d1j+d2jmax*d2j, ...
        d1jmin*d1j+d2jmin*d2j, d1jmax*d1j+d2jmin*d2j];
    
    rectinfo(:,:,i) = [d1j,d2j,[d1jmin;d1jmax],[d2jmin;d2jmax]];
    
end

end

function [rects,rectinfo] = bounding_rects_cheap(convreg,d1s,d2s)
% find tight rectangle for each region given axes

[dim,m,n] = size(convreg);
xreg = reshape(convreg(1,:,:),m,n);
yreg = reshape(convreg(2,:,:),m,n);

rects = zeros(2,4,n);
rectinfo = zeros(2,4,n);

for i = 1:n
    d1i = d1s(:,i);
    d2i = d2s(:,i);
    pts = convreg(:,:,i).';
    d1c = pts*d1i;
    d2c = pts*d2i;
    d1jmax = max(d1c);
    d1jmin = min(d1c);
    d2jmax = max(d2c);
    d2jmin = min(d2c);
    d1j = d1i;
    d2j = d2i;
    
    rects(:,:,i) = [d1jmax*d1j+d2jmax*d2j, d1jmin*d1j+d2jmax*d2j, ...
        d1jmin*d1j+d2jmin*d2j, d1jmax*d1j+d2jmin*d2j];
    
    rectinfo(:,:,i) = [d1j,d2j,[d1jmin;d1jmax],[d2jmin;d2jmax]];
    
end

end
