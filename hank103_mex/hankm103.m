function [h0,h1] = hankm103(z)
%
%     
  asize = size(z);
  tsize = numel(z);
  
  
  h0 = complex(zeros(asize));
  h1 = complex(zeros(asize));

  mex_id_ = 'hank103_wrap(i dcomplex[x], io dcomplex[x], io dcomplex[x], i int64_t[x])';
[h0, h1] = hank103_jgh(mex_id_, z, h0, h1, tsize, tsize, tsize, tsize, 1);

  h0 = reshape(h0,asize);
  h1 = reshape(h1,asize);
end
%
%
%------------------------------------------------------
