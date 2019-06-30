
function [Kpxy,nbr] = proxyfunr(rc,rx,slf,nbr,l,ctr,chnkr,whts,kern,opdims, ...
    pr,ptau,pw,pin)
%PROXYFUNR proxy function utility for kernels defined on chunkers,
% rectangular version
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

pxy = bsxfun(@plus,pr*l,ctr(:));
pw = l*pw;
pw2 = repmat(pw(:).',opdims(1),1); pw2 = pw2(:);

if strcmpi(rc,'c')

    % find unique underlying points corresponding to slf indices

    slfpts = idivide(int64(slf(:)-1),int64(opdims(2)))+1;
    [slfuni,~,islfuni] = unique(slfpts);
    islfuni2 = (islfuni-1)*opdims(2) + mod(slf(:)-1,opdims(2))+1;

    % get matrix-valued entries of kernel for unique points

    rslf = chnkr.r(:,slfuni);
    dslf = chnkr.d(:,slfuni);
    dslf = bsxfun(@rdivide,dslf,sqrt(sum(dslf.^2,1)));

    Kpxy = kern(rslf,pxy,dslf,ptau);

    Kpxy = Kpxy(:,islfuni2);
    Kpxy = bsxfun(@times,Kpxy,whts(slfpts).');

    % targets corresponding to nbr indices

    nbrpts = idivide(int64(nbr(:)-1),int64(opdims(1)))+1;
    dr = (rx(:,nbrpts) - ctr(:))/l;
    nbr = nbr(pin(dr));
else
    % find unique underlying points corresponding to slf indices

    slfpts = idivide(int64(slf(:)-1),int64(opdims(1)))+1;
    [slfuni,~,islfuni] = unique(slfpts);
    islfuni2 = (islfuni-1)*opdims(1) + mod(slf(:)-1,opdims(1))+1;

    % get matrix-valued entries of kernel for unique points

    rslf = rx(:,slfuni);
    dslf = zeros(size(rslf));
    %dslf = bsxfun(@rdivide,dslf,sqrt(sum(dslf.^2,1)));

    Kpxy = kern(pxy,rslf,ptau,dslf);

    Kpxy = Kpxy(islfuni2,:);
    Kpxy = bsxfun(@times,Kpxy,pw2(:).');
    
    % targets corresponding to nbr indices

    nbrpts = idivide(int64(nbr(:)-1),int64(opdims(2)))+1;
    dr = (chnkr.r(:,nbrpts) - ctr(:))/l;
    nbr = nbr(pin(dr));
    
end


end
