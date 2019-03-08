
function [fvals] = curvebymode(t,modes,ctr)
%CURVEBYMODE evaluate the position, and first and second derivatives
% of the position described by r(t) = ctr + (modes(1) + modes(2)*cos(t) +
% modes(3)*sin(t) + ...)(cos(t),sin(t))

  if nargin < 3
    ctr = zeros(2,1);
  end

  x0 = ctr(1); y0 = ctr(2);
  
  eitfac = ones(size(t));
  r = modes(1)*ones(size(t));
  rp = 0.0; rpp = 0.0;
  eit = exp(1i*t);
  for i = 2:2:length(modes)
    eitfac = eitfac.*eit;
    ih = i/2;
    r = r + real(eitfac)*modes(i);
    rp = rp - ih*imag(eitfac)*modes(i);
    rpp = rpp - ih^2*real(eitfac)*modes(i);
    if i < length(modes)
      r = r + imag(eitfac)*modes(i+1);
      rp = rp + ih*real(eitfac)*modes(i+1);
      rpp = rpp - ih^2*imag(eitfac)*modes(i+1);
    end
  end

  xs = x0+r.*real(eit);
  ys = y0+r.*imag(eit);
  dxdr = real(eit);
  dydr = imag(eit);

  dxdth = -r.*imag(eit);
  dxdth2 = -r.*real(eit);
  dydth = r.*real(eit);
  dydth2 = -r.*imag(eit);

  dxdrdth = -imag(eit);
  dydrdth = real(eit);

  dxs = dxdr.*rp+dxdth;
  dys = dydr.*rp+dydth;

  d2xs = dxdr.*rpp+(dxdrdth*2.0d0).*rp+dxdth2;
  d2ys = dydr.*rpp+(dydrdth*2.0d0).*rp+dydth2;
  
  fvals = zeros(length(t),6);


  fvals(:,1) = xs;
  fvals(:,2) = ys;
  fvals(:,3) = dxs;
  fvals(:,4) = dys;
  fvals(:,5) = d2xs;
  fvals(:,6) = d2ys;

end

