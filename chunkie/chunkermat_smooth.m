function mat = chunkermat_smooth(chnkrsrc,chnkrtarg,kern,opdims,glwts)


  targinfo = []; targinfo.r = chnkrtarg.r(:,:); 
  targinfo.d = chnkrtarg.d(:,:); targinfo.d2 = chnkrtarg.d2(:,:);

  srcinfo = []; srcinfo.r = chnkrsrc.r(:,:); 
  srcinfo.d = chnkrsrc.d(:,:); srcinfo.d2 = chnkrsrc.d2(:,:);

  %mat = kern(srcinfo,targinfo);
  hs = chnkrsrc.h;
  dsnrms = sqrt(sum(chnkrsrc.d.^2,1));
  ws = kron(hs(:),glwts(:));
  wts = dsnrms(:).*ws;
  
  %wts = weights(chnkrsrc);
 
  wts2 = repmat(wts(:).',opdims(2),1);
  wts2 = wts2(:);

  mat = bsxfun(@times,kern(srcinfo,targinfo),(wts2(:)).');

end