function flag = flagnear(chnkr,pts,opts)
%FLAGNEAR flag points which are within a chunklength (or given factor of
% a chunklength) of any point on each chunk. on return is a sparse
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
%       opts.fac = factor of chunklength to check distance against (1.0)
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

fac = 1.0;

if nargin < 3
    opts = [];
end

if isfield(opts,'fac')
    fac = opts.fac;
end

npt = chnkr.npt;
nch = chnkr.nch;
k = chnkr.k;
dim = chnkr.dim;

[~,nt] = size(pts);
xx = zeros(2,npt + nt);
xx(:,1:npt) = chnkr.r(:,:);
xx(:,npt+1:npt+nt) = pts;

pt2chnk = repmat(1:nch,k,1);
pt2chnk = pt2chnk(:);
lens = zeros(nch,1);
ws = weights(chnkr);
lens(:) = sum(ws,1); lens = lens*fac;
lmax = max(lens)*2.5;

T = hypoct_uni(xx,lmax);

nnzero = 0;
nn = 3*nt;
isp = zeros(nn,1);
jsp = zeros(nn,1);

for i = 1:length(T.nodes)
    xi = T.nodes(i).xi;
    if nnz(xi <= npt) > 0
        isrc = xi(xi <= npt);
        inbor = [T.nodes(T.nodes(i).nbor).xi T.nodes(i).xi];
        if nnz(inbor > npt) > 0
            
            % get all pairwise distances
            itarg = inbor(inbor > npt)-npt;
            srcs = reshape(chnkr.rstor(:,isrc),dim,1,length(isrc));
            targs = reshape(pts(:,itarg),dim,length(itarg),1);
            dists = reshape(sqrt(sum((bsxfun(@minus,targs,srcs)).^2,1)), ...
                length(itarg),length(isrc));
            
            lenrel = lens(pt2chnk(isrc)); lenrel = (lenrel(:)).';
            distrel = dists < lenrel;
            [it,jt] = find(distrel);
            it = itarg(it);
            jt = pt2chnk(isrc(jt));
            [~,inds] = unique([it(:) jt(:)],'rows');
            it = it(inds); jt = jt(inds);
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
            jsp(nnzero+1:nnzero+nnew) = jt;
            nnzero = nnew + nnzero;
        end
    end
end

isp = isp(1:nnzero); jsp = jsp(1:nnzero);
[~,inds] = unique([isp(:),jsp(:)],'rows');
isp = isp(inds);
jsp = jsp(inds);
vsp = true(length(inds),1);

flag = sparse(isp,jsp,vsp,nt,nch);

end