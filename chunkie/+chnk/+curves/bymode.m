
function [r,d,d2] = bymode(t,modes,ctr,sc)
%CHNK.CURVES.BYMODE evaluate the position, and first and second derivatives
% of the position described by r(t) = ctr + (modes(1) + modes(2)*cos(t) +
% modes(3)*sin(t) + modes(4)*cos(2t) + modes(5)*sin(2t) ...)[cos(t),sin(t)]

  if nargin < 3
    ctr = zeros(2,1);
    sc = [1;1];
  end
  if(nargin < 4)
      sc = [1,1];
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

  xs = x0+r.*real(eit)*sc(1);
  ys = y0+r.*imag(eit)*sc(2);
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
  
  r = [(xs(:)).'; (ys(:)).'];
  d = [(sc(1)*dxs(:)).'; (sc(2)*dys(:)).'];
  d2 = [(sc(1)*d2xs(:)).'; (sc(2)*d2ys(:)).'];  

end

