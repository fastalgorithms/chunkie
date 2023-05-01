% TODO: rewrite this documentation with support for multiple chunkers
function mat = kernbyindex(i,j,chnkrs,whts,kern,opdims_mat,spmat,l2scale)
%% evaluate system matrix by entry index utility function for 
% general kernels, with replacement for specific entries and the
% ability to add a low-rank modification.
%
% opdims_mat:
% We assume that the kernel(i,j) for a given source chunkr(j) and target chunkr(i) is an 
% opdims(1,i,j) x opdims(2,i,j) matrix. 
% 
% Matrix entry indices
% In case with only have 1 chunker object, the matrix entry indices use the following convention: 
% index i corresponds to the (i-1)/opdims(1) + 1 boundary point (integer division implied). 
% Similarly index j corresponds 
% to the (j-1)/opdims(2) + 1 boundary point. 
% For multiple chunker objects, this matrix entry can be generalized
% TODO: elaborate on this
% 
% this function primarily calls the provided kernel routine
% but makes some attempts at efficiency (kernel calls corresponding to
% the same point are not repeated).
%
% input:
%
% i - array of row indices to compute
% j - array of col indices to compute
% chnkr - chunker object describing boundary, 
%         or chunkgraph.echunks, which is an array of chunker objects.
% whts - smooth integration weights on chnkr
% kern - kernel function of the form kern(s,t,stau,ttau) where s and t 
%    are source and target points and stau and ttau are the local unit
%    tangents
% opdims_mat - dimensions of operators, opdims_mat(2,i,j) is the dimension of the 
% spmat - sparse matrix, any entry in output mat corresponding to a 
%    non-zero (non-empty) entry in the matlab built-in sparse 
%    (chnkr.sparse) matrix spmat is overwritten
% see also 
% l2scale = false

% find unique underlying points


nchunkers = length(chnkrs);
irowlocs = zeros(nchunkers+1,1);
icollocs = zeros(nchunkers+1,1);
lchunks = zeros(nchunkers,1);
irowlocs(1) = 1;
icollocs(1) = 1;

for ic=1:nchunkers
   lchunks(ic) = chnkrs(ic).npt;
   icollocs(ic+1) = icollocs(ic) + lchunks(ic)*opdims_mat(2,1,ic);
   irowlocs(ic+1) = irowlocs(ic) + lchunks(ic)*opdims_mat(1,ic,1);
end    

nrows = irowlocs(end)-1;
ncols = icollocs(end)-1;

mat = zeros(length(i),length(j));

if length(i)*length(j) == 0
    return
end

for itrg=1:nchunkers
    chnkr_trg = chnkrs(itrg);
    mat_row_start = irowlocs(itrg);
    mat_row_end   = irowlocs(itrg+1)-1;
    flag_i = i>=mat_row_start & i<=mat_row_end;
    mat_row_ind = i(flag_i);

    if length(mat_row_ind) == 0
        continue
    end

    for isrc=1:nchunkers
        chnkr_src = chnkrs(isrc);
        mat_col_start = icollocs(isrc);
        mat_col_end   = icollocs(isrc+1)-1;
        flag_j = j>=mat_col_start & j<=mat_col_end;
        mat_col_ind   = j(flag_j);

        if length(mat_col_ind) == 0
            continue
        end

        opdims = opdims_mat(:,itrg,isrc);
        ipts = idivide(int64(mat_row_ind(:)-mat_row_start),int64(opdims(1)))+1;
        jpts = idivide(int64(mat_col_ind(:)-mat_col_start),int64(opdims(2)))+1;

        [iuni,~,iiuni] = unique(ipts);
        [juni,~,ijuni] = unique(jpts);

        ri = chnkr_trg.r(:,iuni); rj = chnkr_src.r(:,juni);
        di = chnkr_trg.d(:,iuni); dj = chnkr_src.d(:,juni);
        nj = chnkr_src.n(:,juni); ni = chnkr_trg.n(:,iuni);
        d2i = chnkr_trg.d2(:,iuni); d2j = chnkr_src.d2(:,juni);
        srcinfo = []; srcinfo.r = rj; srcinfo.d = dj; srcinfo.d2 = d2j;
        srcinfo.n = nj;
        targinfo = []; targinfo.r = ri; targinfo.d = di; targinfo.d2 = d2i;
        targinfo.n = ni;
        
        if length(kern) == 1
            matuni = kern(srcinfo,targinfo);
        else
            matuni = kern{itrg,isrc}(srcinfo,targinfo);
        end
        iiuni2 = (iiuni-1)*opdims(1) + mod(mat_row_ind(:)-mat_row_start,opdims(1))+1;
        ijuni2 = (ijuni-1)*opdims(2) + mod(mat_col_ind(:)-mat_col_start,opdims(2))+1;
        mat(flag_i,flag_j) = matuni(iiuni2,ijuni2);
    end
end

% scale columns by weights

if l2scale    
    
    mat = sqrt(whts(i)) .* mat .* sqrt(whts(j).');
    
else
    mat = mat .* whts(j).';
end

% overwrite any entries given as nonzeros in sparse matrix

if nargin > 6
    [isp,jsp,vsp] = find(spmat(i,j));
    linsp = isp + (jsp-1)*length(i(:));
    mat(linsp) = vsp;
end

return;

ipts = idivide(int64(i(:)-1),int64(opdims(1)))+1;
jpts = idivide(int64(j(:)-1),int64(opdims(2)))+1;

[iuni,~,iiuni] = unique(ipts);
[juni,~,ijuni] = unique(jpts);

% matrix-valued entries of kernel for unique points

ri = chnkr.r(:,iuni); rj = chnkr.r(:,juni);
di = chnkr.d(:,iuni); dj = chnkr.d(:,juni);
nj = chnkr.n(:,juni); ni = chnkr.n(:,iuni);
d2i = chnkr.d2(:,iuni); d2j = chnkr.d2(:,juni);
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

whts = whts(:);

% scale columns by weights

if l2scale
    mat = sqrt(whts(ipts(:))) .* mat .* sqrt(whts(jpts(:)).');
else
    mat = mat .* whts(jpts(:)).';
end
% overwrite any entries given as nonzeros in sparse matrix

if nargin > 6
    [isp,jsp,vsp] = find(spmat(i,j));
    linsp = isp + (jsp-1)*length(i(:));
    mat(linsp) = vsp;
end

end
