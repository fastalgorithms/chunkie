function [auxquad]=setup(k,type)
%CHNK.QUADGGQ.SETUP
%
% obtain quadrature nodes and weights for common integral kernel 
% singularities, as well as interpolation matrices for working with these
% nodes. 
%
% quadratures will integrate functions of the form 
%
%     f_0(x) + log(x-x_j)*f_1(x) + 1/(x-x_j)*f_2(x) + 1/(x-x_j)^2*f_3(x)
% [[[---------------------------]-------------------]-------------------]
% type =    'log'                       'pv'                'hs'
%
%     f_{0}(x) + abs|x-x_{j}|*f_{1}(x)
% type = 'removable'
%
% where the f_i are polynomials of order 2*k
%
% input:
%   k - order of Legendre nodes, corresponding to targets x_j 
%   type - the type of singularities to integrate. must be one of 'log', 
%               'pv', or 'hs'
% output:
%   auxquad - struct with entries 
%       auxquad.xs0 - auxquad.xs0{j} are the support nodes for
%       singularities centered at x_j 
%       auxquad.wts0 - auxquad.wts0{j} are the corresponding weights
%       auxquad.ainterps0 - auxquad.ainterps0{j} is an interpolation matrix
%       of size length(auxquad.xs0{j}) x k, interpolating from the order k
%       Legendre nodes to the support nodes.
%       auxquad.xs1 - support nodes for integration on neighbors 
%       auxquad.wts1 - weights for integration on neighbors
%       auxquad.ainterps1 - is an interpolation matrix
%       of size length(auxquad.xs1) x k, interpolating from the order k
%       Legendre nodes to the support nodes for neighbor interactions
%
% NB: the neighbor rule supplied is currently the log type in all cases.

  npolyfac=2;  
  [xs1,wts1,xs0,wts0] = chnk.quadggq.getlogquad(k,npolyfac);
  if strcmpi(type,'pv')
      [xs0,wts0] = chnk.quadggq.gethqsuppquad(k,1);
  elseif strcmpi(type,'hs')
      [xs0,wts0] = chnk.quadggq.gethqsuppquad(k,2);
  elseif strcmpi(type,'removable')
      [xs0,wts0] = chnk.quadggq.getremovablequad(k,1);
  end
  ainterp1 = lege.matrin(k,xs1);

  ainterps0kron = cell(k,1);
  ainterps0 = cell(k,1);
  
  for j = 1:k
    xs0j = xs0{j};
    ainterp0_sm = lege.matrin(k,xs0j);
    ainterps0{j} = ainterp0_sm;
  end
  
  auxquad = [];
  auxquad.xs1 = xs1;
  auxquad.wts1 = wts1;
  auxquad.xs0 = xs0;
  auxquad.wts0 = wts0;
  auxquad.ainterp1 = ainterp1;
  auxquad.ainterps0 = ainterps0;
  
