function [r,d,d2] = chunkerarcparam(s,param_data)
% CHUNKERARCPARAM evaluate arclength parameterization of a chunker
%
% get r, d, d2 at an arbitrary arclength point in a chunker
%
% s is the arclength and should lie in [0, sum(chnkr.wts)]
%
% param_data contains precomputed parameterization data (see
%                           chunkerarcparam_init)
%   param_data.pstrt - arclength coordinate of left endpoints of all
%           panels
%   param_data.plen - arclength of all panels
%   param_data.cr - arclength legendre series of chhnkr.r
%   param_data.cd - arclength legendre series of chhnkr.d
%   param_data.cd2 - arclength legendre series of chhnkr.d2
%   param_data.k - chunker order
%   param_data.dim - chunker dimension
%   param_data.nch - number of chunks in chunkr
%
% see also CHUNKERARCPARAM_INIT

% author: Tristan Goodwill
    
% unpack data
pstrt = param_data.pstrt;
plen = param_data.plen;
cr  = param_data.cr;
cd  = param_data.cd;
cd2 = param_data.cd2;
k = param_data.k;
dim = param_data.dim;
nch = param_data.nch;

% initialize variables
ns = length(s);
r  = zeros(dim, ns);
d  = zeros(dim, ns);
d2 = zeros(dim, ns);

for i = 1:nch
    cri  = cr(:,:,i);
    cdi  = cd(:,:,i);
    cd2i = cd2(:,:,i);

    % determine points in this chunker
    iiin = (s >= pstrt(i)) & (s <= pstrt(i+1));

    sloc = 2*(s(iiin)-pstrt(i))/plen(i) - 1;

    if ~isempty(sloc)
        % evaluate legendre series
        legs = lege.pols(sloc,k-1).';
        r(:,iiin)   = (legs*cri).';
        d(:,iiin)   = (legs*cdi).';
        d2(:,iiin)  = (legs*cd2i).';
    end    
end

end