function [ts,whts] = rts_stab(n)
%LEGE.RTS_STAB computes the nodes and weights for Legendre polynomials
% in O(n^2) time but is generally more stable than LEGE.RTS. 
% A transcription of the original Rokhlin routine
%
% Note: the original Fortran routine supports extended precision, 
% this routine does not
%
%                 Input parameters:
%
%  n - the number of nodes to be returned
%   
%                 Output parameters:
%
%  ts - the n Gaussian nodes on the interval [-1,1]
%  whts - the n Gaussian weights on the interval [-1,1]
% 

% Copyright (C) 2009: Vladimir Rokhlin
% 
% This software is being released under a modified FreeBSD license
%

%
%
%        . . . determine the number of Taylor coefficients 
%              to be used
%


nnewt = 10;
nstop = 3;
stoptol = 1e-12;

% 
%       . . . construct the array of initial approximations
%             to the roots of the n-th legendre polynomial
% 

ifodd = mod(n,2);
%
h=pi/(2*n);

% work on right half of interval, starting guess is Chebyshev
% 
rstart = (n-ifodd)/2+1;
tsr = -cos((2*(rstart:n)-1)*h);
tsr = tsr(:);

%  find roots with Newton

% evaluate at first point O(n), remainder with updates


for kk = 1:length(tsr)
    x1 = tsr(kk);
    
    % conduct newton

    ifstop=0;
    for i=1:nnewt

        [pol,der] = lege.pol(x1,n);
        x1=x1-pol/der;

        if(abs(pol) < stoptol); ifstop=ifstop+1; end
        if(ifstop >= nstop); break; end
    end
    
    tsr(kk)=x1;
end

ts = zeros(n,1);
ts(rstart:n) = tsr(:);
ts(1:(rstart+ifodd-1)) = -flipud(tsr);

if nargout > 1
    [~,~,totr] = lege.polsum(tsr,n);
    whts = zeros(n,1);
    whts(rstart:n) = 1./totr(:);
    whts(1:(rstart+ifodd-1)) = flipud(1./totr(:));
end
