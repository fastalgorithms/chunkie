function mat = chunkermat_smooth(chnkrsrc,chnkrtarg,kern)


  targinfo = []; targinfo.r = chnkrtarg.r(:,:); 
  targinfo.d = chnkrtarg.d(:,:); targinfo.d2 = chnkrtarg.d2(:,:);
  nt = chnkrtarg.k*chnkrtarg.nch;

  srcinfo = []; srcinfo.r = chnkrsrc.r(:,:); 
  srcinfo.d = chnkrsrc.d(:,:); srcinfo.d2 = chnkrsrc.d2(:,:);

  mat = kern(srcinfo,targinfo);

  wts = weights(chnkrsrc);
  wts2 = repmat( (wts(:)).', nt, 1);

  mat = mat.*wts2;

end