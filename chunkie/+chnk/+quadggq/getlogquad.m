function [xs1,ws1,xs0,ws0] = getlogquad(k,npolyfac)
%GETLOGQUAD returns a GGQ quadrature for smooth+log*smooth
% of the requested order (if available)
%
% input
% k - order of nodes on panel
% npolyfac - number of polynomials to integrate is k*npolyfac
%             (must be either 2 or 3)

  if (nargin < 2)
    npolyfac = 2;
  end
  
  if k <= 16
      kmax = 16;
  elseif 16 < k && k <= 20
      kmax = 20;
  elseif 20 < k && k <= 24
      kmax = 24;
  elseif 24 < k && k <= 30
      kmax = 30;
  elseif 30 < k && k <= 40
      kmax = 40;
  elseif 40 < k && k <= 60
      kmax = 60;
  else
      error('quadggq: selected order has no corresponding near rule');
  end
  
  [xs1,ws1] = chnk.quadggq.ggqnear(kmax);

  switch k
    case 1
      [xs0,ws0] = chnk.quadggq.ggqself_nnode001_npoly002();
    case 2
      [xs0,ws0] = chnk.quadggq.ggqself_nnode002_npoly004();
    case 3
      [xs0,ws0] = chnk.quadggq.ggqself_nnode003_npoly006();
    case 4
      [xs0,ws0] = chnk.quadggq.ggqself_nnode004_npoly008();
    case 5
      [xs0,ws0] = chnk.quadggq.ggqself_nnode005_npoly010();
    case 6
      [xs0,ws0] = chnk.quadggq.ggqself_nnode006_npoly012();
    case 7
      [xs0,ws0] = chnk.quadggq.ggqself_nnode007_npoly014();
    case 8
      [xs0,ws0] = chnk.quadggq.ggqself_nnode008_npoly016();
    case 9
      [xs0,ws0] = chnk.quadggq.ggqself_nnode009_npoly018();
    case 10
      [xs0,ws0] = chnk.quadggq.ggqself_nnode010_npoly020();
    case 11
      [xs0,ws0] = chnk.quadggq.ggqself_nnode011_npoly022();
    case 12
      [xs0,ws0] = chnk.quadggq.ggqself_nnode012_npoly024();
    case 13
      [xs0,ws0] = chnk.quadggq.ggqself_nnode013_npoly026();
    case 14
      [xs0,ws0] = chnk.quadggq.ggqself_nnode014_npoly028();
    case 15
      [xs0,ws0] = chnk.quadggq.ggqself_nnode015_npoly030();
    case 16
      [xs0,ws0] = chnk.quadggq.ggqself_nnode016_npoly032();
    case 20
      [xs0,ws0] = chnk.quadggq.ggqself_nnode020_npoly040();
    case 24
      [xs0,ws0] = chnk.quadggq.ggqself_nnode024_npoly048();
    case 28
      [xs0,ws0] = chnk.quadggq.ggqself_nnode028_npoly056();
    case 32
      [xs0,ws0] = chnk.quadggq.ggqself_nnode032_npoly064();
    case 36
      [xs0,ws0] = chnk.quadggq.ggqself_nnode036_npoly072();
    %case 40
    %  [xs0,ws0] = chnk.quadggq.ggqself_nnode040_npoly080();
    otherwise
      error("getlogquad: order not available\n");
  end	
end
