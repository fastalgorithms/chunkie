function [ts,whts] = rts(n)
%LEGE.RTS computes the nodes and weights for Legendre polynomials
% in O(n) time. A transcription of the original Rokhlin routine
%
%        This subroutine constructs the Gaussian quadrature
%        or order n. Its claim to fame is the fact that the
%        cost of the calculation is proportional to n; in
%        practice, with n=10 000 the calculation is more or
%        less instantaneous
%
% Note: the original Fortran routine supports extended precision, 
% this routine does not
%
%                 Input parameters:
%
%  itype - the type of calculation desired: 
%     itype=1 will cause both the roots and the weights to be returned
%     itype=0 will cause only the roots to be returned
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

k = 30; % double precision

% extended precision not supported
%d=1.0;
%d2=d+1.0e-24;
%if (d2 ~= d); k=54; end

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
whtsr = zeros(length(tsr),1);

%  find roots with Newton

% evaluate at first point O(n), remainder with updates

x0 = 0.0;
[pol,der] = lege.pol(x0,n); 

polpre = pol;
derpre = der;

for kk = 1:length(tsr)
    % for odd, zero is correct
    if( and(ifodd == 1, kk == 1) ) 
        tsr(kk)=0.0;
        whtsr(kk)=der;
        if (n==0); return; end
        polpre=pol;
        derpre=der;
        continue
    end

    x1 = tsr(kk);
    
    % conduct newton

    ifstop=0;
    for i=1:nnewt

        h=x1-x0;
        [pol,der] = lege.tayl(polpre,derpre,x0,h,n,k);
        x1=x1-pol/der;

        if(abs(pol) < stoptol); ifstop=ifstop+1; end
        if(ifstop >= nstop); break; end
    end

    tsr(kk)=x1;
    whtsr(kk)=der;
    x0=x1;
    polpre=pol;
    derpre=der;
end

ts = zeros(n,1);
ts(rstart:n) = tsr(:);
ts(1:(rstart+ifodd-1)) = -flipud(tsr);

if nargout > 1
    whts(rstart:n) = whtsr(:);
    whts(1:(rstart+ifodd-1)) = flipud(whtsr);
    whts = 2.0./(1-ts(:).^2)./whts(:).^2;
end
  

end


