
function [proxy,pnorm,pw] = proxy_square_pts(p)

if nargin < 1
    p = 64;
end

po4 = p/4;

pts = -1.5+3*(0:(po4-1))*1.0/po4;
one = ones(size(pts));

proxy = [pts, one*1.5, -pts, -1.5*one;
    -1.5*one, pts, 1.5*one, -pts];

pnorm = [one*0, one, one*0, -one;
    -one, one*0, one, one*0];

pw = ones(p,1)*(4.0*3.0/p);

end