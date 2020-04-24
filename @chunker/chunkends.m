function [rend,tauend] = chunkends(chnkr,ich)
%CHUNKENDS evaluate the endpoints of chunks in chunker object
%
% Syntax: [rend,tauend] = ends(chnkr,ich)
%
% Input:
%   chnkr - chunker object
%
% Optional input:
%   ich - list of indices (if not specified, computes all ends)
%
% Output:
%   rend - (dim,2,:) array, positions of end points for requested chuns
%   tauend - (dim,2,:) array, unit tangent vectors at end points
%
% Examples:
%   rend = ends(chnkr);
%   ich = [1 3 9];
%   [rend,tauend] = ends(chnkr,ich);
%

% author: Travis Askham (askhamwhat@gmail.com)

[~,~,u] = lege.exps(chnkr.k);
pm1 = lege.pols(-1,chnkr.k-1); pm1 = pm1.';
p1 = lege.pols(1,chnkr.k-1); p1 = p1.';

pends = [pm1; p1];

if nargin < 2
    nch1 = chnkr.nch;
    dim = chnkr.dim;
    ich = 1:nch1;
else
    nch1 = numel(ich);
end

dim = chnkr.dim;
rend = zeros(dim,2,nch1);
tauend = zeros(dim,2,nch1);

for i = 1:nch1
    ii = ich(i);
    ri = chnkr.r(:,:,ii);
    di = chnkr.d(:,:,ii);    
    cri = u*(ri.');
    cdi = u*(di.');
    rendi = pends*cri; rendi = rendi.';
    dendi = pends*cdi; dendi = dendi.';
    
    dsendi = sqrt(sum(dendi.^2,1));
    tauendi = bsxfun(@rdivide,dendi,dsendi);
    
    rend(:,:,i) = rendi;
    tauend(:,:,i) = tauendi;
end
    
