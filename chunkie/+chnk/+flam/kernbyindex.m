function mat = kernbyindex(i,j,chnkr,whts,kern,opdims,spmat)
%% evaluate system matrix by entry index utility function for 
% general kernels, with replacement for specific entries and the
% ability to add a low-rank modification.
%
% We assume that the kernel for a given source and target is an 
% opdims(1) x opdims(2) matrix. Matrix entry indices use the 
% following convention: index i corresponds to the (i-1)/opdims(1) + 1 
% boundary point (integer division implied). Similarly index j corresponds 
% to the (j-1)/opdims(2) + 1 boundary point. 
%
% this function primarily calls the provided kernel routine
% but makes some attempts at efficiency (kernel calls corresponding to
% the same point are not repeated).
%
% input:
%
% i - array of row indices to compute
% j - array of col indices to compute
% chnkr - chunker object describing boundary
% whts - smooth integration weights on chnkr
% kern - kernel function of the form kern(s,t,stau,ttau) where s and t 
%    are source and target points and stau and ttau are the local unit
%    tangents
% opdims - dimensions of operator (opdims(1) dim of output, opdims(2) dim
%    of input)
% spmat - sparse matrix, any entry in output mat corresponding to a 
%    non-zero (non-empty) entry in the matlab built-in sparse 
%    (chnkr.sparse) matrix spmat is overwritten
% see also 

% find unique underlying points

ipts = idivide(int64(i(:)-1),int64(opdims(1)))+1;
jpts = idivide(int64(j(:)-1),int64(opdims(2)))+1;

[iuni,~,iiuni] = unique(ipts);
[juni,~,ijuni] = unique(jpts);

% matrix-valued entries of kernel for unique points

ri = chnkr.r(:,iuni); rj = chnkr.r(:,juni);
di = chnkr.d(:,iuni); dj = chnkr.d(:,juni);
nj = chnkr.n(:,juni); ni = chnkr.n(:,iuni);
d2i = chnkr.d(:,iuni); d2j = chnkr.d2(:,juni);
srcinfo = []; srcinfo.r = rj; srcinfo.d = dj; srcinfo.d2 = d2j;
srcinfo.n = nj;
targinfo = []; targinfo.r = ri; targinfo.d = di; targinfo.d2 = d2i;
targinfo.n = ni;
%di = bsxfun(@rdivide,di,sqrt(sum(di.^2,1)));
%dj = bsxfun(@rdivide,dj,sqrt(sum(dj.^2,1)));

matuni = kern(srcinfo,targinfo);

% relevant rows and columns in matuni

iiuni2 = (iiuni-1)*opdims(1) + mod(i(:)-1,opdims(1))+1;
ijuni2 = (ijuni-1)*opdims(2) + mod(j(:)-1,opdims(2))+1;

mat = matuni(iiuni2,ijuni2);

% scale columns by weights

wj = whts(jpts(:));
mat = bsxfun(@times,mat,wj.');

% overwrite any entries given as nonzeros in sparse matrix

if nargin > 6
    [isp,jsp,vsp] = find(spmat(i,j));
    linsp = isp + (jsp-1)*length(i(:));
    mat(linsp) = vsp;
end

end
