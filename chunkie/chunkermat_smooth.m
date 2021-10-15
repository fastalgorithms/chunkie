function mat = chunkermat_smooth(chnkrsrc,chnkrtarg,kern,opdims,glwts,slist,tlist)

  targinfo = []; targinfo.r = chnkrtarg.r(:,:); targinfo.n = chnkrtarg.n(:,:);
  targinfo.d = chnkrtarg.d(:,:); targinfo.d2 = chnkrtarg.d2(:,:);

  srcinfo = []; srcinfo.r = chnkrsrc.r(:,:); srcinfo.n = chnkrsrc.n(:,:);
  srcinfo.d = chnkrsrc.d(:,:); srcinfo.d2 = chnkrsrc.d2(:,:);

  if isempty(slist) && isempty(tlist)

  else
    k = chnkrsrc.k;
    
    is = 1:chnkrsrc.nch;
    it = 1:chnkrtarg.nch;
    
    isg = setdiff(is, slist);
    isb = slist';
    
    itg = setdiff(it, tlist);
    itb = tlist';
    
    sg = []; sg.r = chnkrsrc.r(:,:,isg); sg.n = chnkrsrc.n(:,:,isg);
    sb = []; sb.r = chnkrsrc.r(:,:,isb); sb.n = chnkrsrc.n(:,:,isb);
    
    tg = []; tg.r = chnkrtarg.r(:,:,itg); tg.n = chnkrtarg.n(:,:,itg);
    tb = []; tb.r = chnkrtarg.r(:,:,itb); tb.n = chnkrtarg.n(:,:,itb);
    
    mat1 = kern(sg,tg);
    mat2 = kern(sb,tg);
    mat3 = kern(sg,tb);
    
    ns = opdims(2)*chnkrsrc.nch*k;
    nt = opdims(1)*chnkrtarg.nch*k;
    
    mat = zeros(nt,ns);
    
    sn = opdims(2)*k; tn = opdims(1)*k;
    is0 = 1:sn; it0 = 1:tn;
    sgid = [];sbid=[];tgid=[];tbid=[];
    for i=isg
      sgid = [sgid (i-1)*sn+is0];
    end
    for i=isb
      sbid = [sbid (i-1)*sn+is0];
    end
    for i=itg
      tgid = [tgid (i-1)*tn+it0];
    end
    for i=itb
      tbid = [tbid (i-1)*tn+it0];
    end   
    
    mat(tgid,sgid) = mat1;
    mat(tgid,sbid) = mat2;
    mat(tbid,sgid) = mat3;
  end
    
    
  hs = chnkrsrc.h;
  dsnrms = sqrt(sum(chnkrsrc.d.^2,1));
  ws = kron(hs(:),glwts(:));
  wts = dsnrms(:).*ws;
  
  %wts = weights(chnkrsrc);
 
  wts2 = repmat(wts(:).',opdims(2),1);
  wts2 = wts2(:);
  
  if isempty(slist) && isempty(tlist)
    mat = bsxfun(@times,kern(srcinfo,targinfo),(wts2(:)).');
  else
    mat = bsxfun(@times,mat,(wts2(:)).');
  end
  %max(mat,[],'all')
end