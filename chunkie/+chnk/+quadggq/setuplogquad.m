function [logquad]=setuplogquad(k,opdims)
% obtain quadrature nodes and weights, and interpolation matrices for
% logarithmic and nearly logarithmic singularities

msg = "chnk.quadggq.setuplogquad to be deprecated. " + ...
    "use chnk.quadggq.setup instead";
warning(msg);

  npolyfac=2;  
  [xs1,wts1,xs0,wts0] = chnk.quadggq.getlogquad(k,npolyfac);
  ainterp1 = lege.matrin(k,xs1);
  temp = eye(opdims(2));
  ainterp1kron = kron(ainterp1,temp);

  nquad0 = size(xs0,1);

  ainterps0kron = cell(k,1);
  ainterps0 = cell(k,1);
  
  for j = 1:k
    xs0j = xs0{j};
    ainterp0_sm = lege.matrin(k,xs0j);
    ainterps0{j} = ainterp0_sm;
    ainterps0kron{j} = kron(ainterp0_sm,temp);
  end
  
  logquad = [];
  logquad.xs1 = xs1;
  logquad.wts1 = wts1;
  logquad.xs0 = xs0;
  logquad.wts0 = wts0;
  logquad.ainterp1 = ainterp1;
  logquad.ainterp1kron = ainterp1kron;
  logquad.ainterps0 = ainterps0;
  logquad.ainterps0kron = ainterps0kron;
  
