function in = chunkerinteriorflam(chnkr,pts)
%CHUNKERINTERIORFLAM 

kernd = @(s,t,sn,tn) chnk.lap2d.kern(s,t,sn,tn,'d');
dens1 = ones(chnkr.k,chnkr.nch);
wts = weights(chnkr);

opdims = [1 1];

xflam = chnkr.r(:,:);
matfun = @(i,j) kernbyindexr(i,j,pts,chnkr,wts,kernd,opdims);
[pr,ptau,pw,pin] = proxy_square_pts();

pxyfun = @(rc,rx,cx,slf,nbr,l,ctr) proxyfunr(rc,rx,slf,nbr,l, ...
    ctr,chnkr,wts,kernd,opdims,pr,ptau,pw,pin);
F = ifmm(matfun,pts,xflam,200,1e-14,pxyfun);
vals1 = ifmm_mv(F,dens1(:),matfun);

in = abs(vals1+1) < abs(vals1);

% find nearest neighbors at certain level of refinement (here chosen
% uniformly, for simplicity)

nt = size(pts,2);
xflam = zeros(2,chnkr.npt + nt);
xflam(:,1:chnkr.npt) = chnkr.r(:,:);
xflam(:,chnkr.npt+1:chnkr.npt+nt) = pts;

pt2chnk = repmat(1:chnkr.nch,chnkr.k,1);

chunklens = zeros(chnkr.nch,1);
ws = weights(chnkr);
chunklens(:) = sum(ws,1);
lmax = max(chunklens)*3/chnkr.k;

T = hypoct_uni(xflam,lmax);
targ_dists = Inf(nt,1);
itarg_dists = (1:nt).';

rnorm = normals(chnkr);
normdist_flag = false(size(targ_dists));

for i = 1:length(T.nodes)
    xi = T.nodes(i).xi;
    if nnz(xi <= chnkr.npt) > 0
        isrc = xi(xi <= chnkr.npt);
        inbor = [T.nodes(T.nodes(i).nbor).xi T.nodes(i).xi];
        if nnz(inbor > chnkr.npt) > 0
            itarg = inbor(inbor > chnkr.npt)-chnkr.npt;
            srcs = reshape(chnkr.r(:,isrc),chnkr.dim,1,length(isrc));
            targs = reshape(pts(:,itarg),chnkr.dim,length(itarg),1);
            dists = reshape(sqrt(sum((bsxfun(@minus,targs,srcs)).^2,1)), ...
                length(itarg),length(isrc));
            [mindists,inds] = min(dists,[],2);
            isrc2 = isrc(inds);
            
            ifnew = mindists < targ_dists(itarg);
            targ_dists(itarg(ifnew)) = mindists(ifnew);
            itarg_dists(itarg(ifnew)) = isrc2(ifnew);
            
            rdiff = chnkr.r(:,isrc2(ifnew)) - pts(:,itarg(ifnew));
            rnorms = rnorm(:,isrc2(ifnew));
            lens = chunklens(pt2chnk(isrc2(ifnew)));
            
            rdrn = sum(rdiff.*rnorms,1); rdrn = rdrn(:);
            lens = lens(:);
            in(itarg(ifnew)) = rdrn > 0;
            normdist_flag(itarg(ifnew)) = abs(rdrn) < 0.1*lens;            
        end
    end
end

%

iflagged = find(normdist_flag);
[tt,~,u] = lege.exps(chnkr.k);
for i = 1:length(iflagged)
    ii = iflagged(i);
    ich = pt2chnk(itarg_dists(ii));
    ich = [ich; chnkr.adj(:,ich)];
    [rn,dn] = nearest(chnkr,pts(:,ii),ich,tt,u);
    in(ii) = (rn(:)-pts(:,ii)).'*[dn(2);-dn(1)] > 0;
end