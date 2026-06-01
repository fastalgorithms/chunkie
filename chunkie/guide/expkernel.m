function mat = expkernel(s,t,k)
%EXPKERNEL this implements the kernel evaluator function for the
% kernel exp(ik|x-y|) in the chunkIE kernel format as part of the
% chunkIE guide   
%
% Syntax:
%    mat = expkernel(s,t,k) 
%   

  xs = s.r(1,:);
  ys = s.r(2,:);
  xt = t.r(1,:);
  yt = t.r(2,:);

  mat = exp(1i*k*sqrt((xt.'-xs).^2+(yt.'-ys).^2));
