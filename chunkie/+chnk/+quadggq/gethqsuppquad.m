function [xs0,ws0] = gethqsuppquad(k,itype)
%GETLOGQUAD returns a GGQ quadrature for smooth+log*smooth
% of the requested order (if available)
%
% input
% k - order of nodes on panel
% itype - 1 gives a rule for log and 1/(x-xk), 2 gives a rule
%           for log, 1/(x-xk), 1/(x-xk)^2   

  if nargin < 2
    itype = 2;
  end

  if itype == 1
    switch k
      case 1
	[xs0,ws0] = chnk.quadggq.hsupp_nnode001_npoly002();
      case 2
	[xs0,ws0] = chnk.quadggq.hsupp_nnode002_npoly004();
      case 3
	[xs0,ws0] = chnk.quadggq.hsupp_nnode003_npoly006();
      case 4
	[xs0,ws0] = chnk.quadggq.hsupp_nnode004_npoly008();
      case 5
	[xs0,ws0] = chnk.quadggq.hsupp_nnode005_npoly010();
      case 6
	[xs0,ws0] = chnk.quadggq.hsupp_nnode006_npoly012();
      case 8
	[xs0,ws0] = chnk.quadggq.hsupp_nnode008_npoly016();
      case 10
	[xs0,ws0] = chnk.quadggq.hsupp_nnode010_npoly020();
      case 12
	[xs0,ws0] = chnk.quadggq.hsupp_nnode012_npoly024();
      case 16
	[xs0,ws0] = chnk.quadggq.hsupp_nnode016_npoly032();
      case 20
	[xs0,ws0] = chnk.quadggq.hsupp_nnode020_npoly040();
      otherwise
	error("gethqsuppquad: order not available\n");
    end
  end
  if itype == 2
    switch k
      case 1
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode001_npoly002();
      case 2
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode002_npoly004();
      case 3
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode003_npoly006();
      case 4
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode004_npoly008();
      case 5
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode005_npoly010();
      case 6
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode006_npoly012();
      case 8
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode008_npoly016();
      case 10
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode010_npoly020();
      case 12
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode012_npoly024();
      case 16
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode016_npoly032();
      case 20
	[xs0,ws0] = chnk.quadggq.hqsupp_nnode020_npoly040();
      otherwise
	error("gethqsuppquad: order not available\n");
    end
  end
end
