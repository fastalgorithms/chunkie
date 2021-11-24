function [Kpxy,nbr] = proxyfun_ell2_multi(slf_in,nbr_in,l,ctr,chnkr,cwhts,...
    kern,opdims,pr,ptau,pw,pin,ifaddtrans,ks,opts_perm)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%PROXYFUN proxy function utility for kernels defined on chunkers
%
% We assume that the kernel for a given source and target is an 
% opdims(1) x opdims(2) matrix. Matrix entry indices use the 
% following convention: Similarly index j corresponds 
% to the (j-1)/opdims(2) + 1 boundary point for slf points and 
% index i corresponds to the (i-1)/opdims(1) + 1 boundary point
% for the nbr points (from the point of view of this routine slf 
% points will be sources where nbr points are potential targets)
%
% this function primarily calls the provided kernel routine
% but makes some attempts at efficiency (kernel calls corresponding to
% the same point are not repeated).
%
% ~ octtree and index related inputs ~
% slf - set of relevant indices for self box (not necessarily equivalent
%       to the appropriate indices within the chunker object)
% nbr - set of indices in neighbor boxes (to check for inclusion inside 
%        proxy surface)
% l - box side length at current level
% ctr - center of current box
%
% ~ other inputs ~
% chnkr - chunker object 
% whts - smooth weights on chunker object
% kern - should be a kernel function of the form 
%             submat = kern(src, targ, srctau, targtau)
% opdims - dimensions of operator, opdims(1) output dim, opdims(2) input
%            dim of operator (distinct from spatial dim)
% pr - a (dim,_) array of proxy points (dim is spatial dim here)
% ptau - a (dim,_) array of proxy unit tangents
% pw - is a set of smooth integration weights for the proxy points
% pin - function handle, takes in rescaled and recentered points in the 
%       nbr array and determines if they're within the proxy surface
%

% scaled proxy points and weights (no scaling necessary on tangents)

lmax = max(l);
pxy = bsxfun(@plus,pr*lmax,ctr(:));
pw = lmax*pw;
pw2 = repmat(pw(:).',opdims(1),1); pw2 = pw2(:);

% find unique underlying points corresponding to slf indices
slf = opts_perm.iperm(slf_in);
slfpts = idivide(int64(slf(:)-1),int64(opdims(2)))+1;
[slfuni,~,islfuni] = unique(slfpts);
islfuni2 = (islfuni-1)*opdims(2) + mod(slf(:)-1,opdims(2))+1;

% get matrix-valued entries of kernel for unique points

srcinfo = []; srcinfo.r = chnkr.r(:,slfpts);
srcinfo.n = chnkr.n(:,slfpts);
srcinfo.d = chnkr.d(:,slfpts); srcinfo.d2 = chnkr.d2(:,slfpts);
targinfo = []; targinfo.r = pxy; targinfo.d = ptau; 
targinfo.d2 = []; targinfo.n = ptau;
pw_lng = [pw;pw;pw];
pw_lng = [pw_lng;pw_lng];
pw_lng = sqrt(pw_lng);

Kpxy = [];

for i=1:numel(ks)
zk = ks(i);

Kpxy_nt = kern(zk,srcinfo,targinfo);
%Kpxy_t = Kpxy_t(:,islfuni2,:);
Kpxy_nt = permute(Kpxy_nt,[1,3,2]);
nc = size(Kpxy_nt,3);
nr = numel(Kpxy_nt)/nc;
Kpxy_nt = reshape(Kpxy_nt,[nr,nc]);
Kpxy_nt = bsxfun(@times,Kpxy_nt,sqrt(cwhts(slfpts)).');
Kpxy_nt = bsxfun(@times,Kpxy_nt,pw_lng);

if ifaddtrans
    Kpxy_t = kern(zk,targinfo,srcinfo);
    Kpxy_t = reshape(Kpxy_t,[nc,nr]);
    %Kpxy_t = Kpxy2(islfuni2,:);
    Kpxy_t = bsxfun(@times,Kpxy_t,pw_lng(:).');
    Kpxy_t = bsxfun(@times,Kpxy_t,sqrt(cwhts(slfpts)));
    
end

Kpxy = [Kpxy_nt; Kpxy_t.'];
end
% points of geometry corresponding to nbr indices

nbr = opts_perm.iperm(nbr_in);
nbrpts = idivide(int64(nbr(:)-1),int64(opdims(1)))+1;
dr = (chnkr.r(:,nbrpts) - ctr(:))/lmax;
nbr = nbr_in(pin(dr));

end
