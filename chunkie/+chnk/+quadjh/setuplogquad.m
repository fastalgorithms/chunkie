function [LogC]=setuplogquad(hs,k,isclose,opdims)
% obtain the correction matrix for logarithmic singularities in
% kernel-split quadrature
[xs,ws] = lege.exps(k);
M1 = chnk.quadjh.LogCinit(hs,xs,ws,isclose);
LogC = kron(M1,ones(opdims));

% logquad = [];
% logquad.LogC = LogC;