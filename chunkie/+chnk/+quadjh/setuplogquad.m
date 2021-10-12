function [LogC]=setuplogquad(hs,k,isclosed,opdims)
% obtain the correction matrix for logarithmic singularities in
% kernel-split quadrature
[xs,ws] = lege.exps(k);
M1 = chnk.quadjh.LogCinit(hs,xs,ws,isclosed);
LogC = kron(M1,ones(opdims));

% logquad = [];
% logquad.LogC = LogC;