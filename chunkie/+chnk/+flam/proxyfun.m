% TODO: change the documentation as multipple chunkers are adopted... 

function [Kpxy,nbr] = proxyfun(slf,nbr,l,ctr,chnkrs,kern,opdims_mat, ...
    pr,ptau,pw,pin,ifaddtrans,l2scale)
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
% ~ tree and index related inputs ~
% slf - set of relevant indices for self box (not necessarily equivalent
%       to the appropriate indices within the chunker object)
% nbr - set of indices in neighbor boxes (to check for inclusion inside 
%        proxy surface)
% l - box side length at current level
% ctr - center of current box
%
% ~ other inputs ~
% chnkrs - chunker object or an array of chunker objects
% kern - should be a kernel function of the form 
%             submat = kern(src, targ, srctau, targtau)
%        or multiple kernel functions for the multiple chunker case
% opdims_mat - dimensions of operator, opdims(1) output dim, opdims(2) input
%            dim of operator (distinct from spatial dim)
%            size = (2,nchunker,nchunker)
% pr - a (dim,_) array of proxy points (dim is spatial dim here)
% ptau - a (dim,_) array of proxy unit tangents
% pw - is a set of smooth integration weights for the proxy points
% pin - function handle, takes in rescaled and recentered points in the 
%       nbr array and determines if they're within the proxy surface
% l2scale - boolean type that determines if we should 
%    rescale the matrix by l2scale. the default value is false. 
% 

if ~isa(kern,'kernel')
    try 
        kern = kernel(kern);
    catch
        error('PROXYFUN: sixth input kern not of supported type');
    end
end


% scaled proxy points and weights (no scaling necessary on tangents)

lmax = max(l);
pxy = pr*lmax + ctr(:);
pw = lmax*pw;

npxy = size(pxy,2);

nchunkers = length(chnkrs);
irowlocs = zeros(nchunkers+1,1);
irowlocs_trans = zeros(nchunkers+1,1);
icollocs = zeros(nchunkers+1,1);
inbrlocs = zeros(nchunkers+1,1);
lchunks = zeros(nchunkers,1);
irowlocs(1) = 1;
irowlocs_trans(1) = 1;
icollocs(1) = 1;
inbrlocs(1) = 1;

for i=1:nchunkers
   lchunks(i) = chnkrs(i).npt;
   icollocs(i+1) = icollocs(i) + lchunks(i)*opdims_mat(2,1,i);
   irowlocs(i+1) = irowlocs(i) + npxy*opdims_mat(1,i,1);
   irowlocs_trans(i+1) = irowlocs_trans(i) + npxy*opdims_mat(2,1,i);
   inbrlocs(i+1) = inbrlocs(i) + lchunks(i)*opdims_mat(1,i,1);
end

nrows = irowlocs(end)-1;
nrows_trans = irowlocs_trans(end)-1;
ncols = icollocs(end)-1;


% % points of geometry corresponding to nbr indices

new_nbr = [];

for i=1:nchunkers
    chnkr = chnkrs(i);
    opdim = opdims_mat(1,i,1);
    fnbr = nbr>=inbrlocs(i) & nbr< inbrlocs(i+1);
    inbr = nbr(fnbr);
    if isempty(inbr)
        continue
    end
    nbrpts = idivide(int64(inbr(:)-inbrlocs(i)),int64(opdim))+1;
    dr = (chnkr.r(:,nbrpts) - ctr(:))/lmax;
    new_nbr = [new_nbr, inbr(pin(dr))];
end

nbr = new_nbr;


% % get matrix-valued entries of kernel for unique points

Kpxy = zeros(nrows,length(slf));

if ifaddtrans
    Kpxy2 = zeros(length(slf),nrows_trans);
end

targinfo = []; targinfo.r = pxy; targinfo.d = ptau; 
targinfo.d2 = []; targinfo.n = [-ptau(2,:); ptau(1,:)] ./ sqrt(sum(ptau.^2,1));


for j=1:nchunkers

    chnkr = chnkrs(j);
    mat_col_start = icollocs(j);
    mat_col_end   = icollocs(j+1)-1;
    f_col = slf>=mat_col_start & slf<=mat_col_end;
    mat_col_ind   = slf(f_col);
    wts = chnkr.wts; wts = wts(:);

% Note that opdims(2) must be same for all rows and a fixed column
    opdims = opdims_mat(:,1,j);

    jpts = idivide(int64(mat_col_ind(:)-mat_col_start),int64(opdims(2)))+1;
    [juni,~,ijuni] = unique(jpts);

    ijuni2 = (ijuni-1)*opdims(2) + mod(mat_col_ind(:)-mat_col_start,opdims(2))+1;

    r = chnkr.r(:,juni);
    d = chnkr.d(:,juni);
    n = chnkr.n(:,juni);
    d2 = chnkr.d2(:,juni);
    srcinfo = []; srcinfo.r = r; srcinfo.d = d; srcinfo.n = n; srcinfo.d2 = d2;
    
    if isempty(mat_col_ind)
        continue
    end

    for i=1:nchunkers

        mat_row_start = irowlocs(i);
        mat_row_end = irowlocs(i+1)-1;
        mat_row_ind = mat_row_start:mat_row_end;

        mat_row_trans_start = irowlocs_trans(i);
        mat_row_trans_end = irowlocs_trans(i+1)-1;
        mat_row_trans_ind = mat_row_trans_start:mat_row_trans_end;
        opdims = opdims_mat(:,i,j);
        opdims_trans = opdims_mat(:,j,i);

        if length(kern) == 1
            matuni = kern.eval(srcinfo,targinfo);
        else
            matuni = kern(i,j).eval(srcinfo,targinfo);
        end

        wsrc = wts(juni);
        wsrc = repmat( (wsrc(:)).', opdims(2), 1);
        wsrc = ( wsrc(:) ).';

        wtarg = pw(:);
        wtarg = repmat(wtarg(:).',opdims(1),1);
        wtarg = wtarg(:);
        
        
        % scale matrix by weights
        if(l2scale)
            matuni = sqrt(wtarg) .* matuni .* sqrt(wsrc);
        else
            matuni = matuni .* wsrc;
        end
        


        Kpxy(mat_row_ind,f_col) = matuni(:,ijuni2);

        if ifaddtrans
            if length(kern) == 1
                matuni2 = kern.eval(targinfo,srcinfo);
            else
                matuni2 = kern(j,i).eval(targinfo,srcinfo);
            end

            wsrc = pw(:);
            wsrc = repmat( (wsrc(:)).', opdims_trans(2), 1);
            wsrc = ( wsrc(:) ).';

            wtarg = wts(juni);
            wtarg = repmat(wtarg(:).',opdims_trans(1),1);
            wtarg = wtarg(:);

            % scale matrix by weights
            if(l2scale)
                matuni2 = sqrt(wtarg) .* matuni2 .* sqrt(wsrc);
            else
                matuni2 = matuni2 .* wsrc;
            end
        

            Kpxy2(f_col,mat_row_trans_ind) = matuni2(ijuni2,:);
        end
    end
end


% combining the two matrices
if ifaddtrans
    Kpxy = [Kpxy; Kpxy2.'];
end

end
