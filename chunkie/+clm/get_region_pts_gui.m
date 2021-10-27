function [x] = get_region_pts_gui(chnkr,clmparams,ireg)
    nch = 0;
    k = chnkr(1).k;
    iregstart = 1;
    iregend = length(clmparams.clist{ireg});
    iextra = 1;
    istart = 1;
    if(ireg == 1 || ireg == 2)
        iregstart = 2;
        iregend = iregend - 1;
        iextra = 3;
        istart = 2;
    end
    for i=iregstart:iregend
        nch = nch + chnkr(abs(clmparams.clist{ireg}(i))).nch;
    end
    
    
    x = zeros(2,nch*k+iextra);
    
    for i=iregstart:iregend
        rtmp = chnkr(abs(clmparams.clist{ireg}(i))).r;
        nch0 = chnkr(abs(clmparams.clist{ireg}(i))).nch;
        rtmp = reshape(rtmp,[2,k*nch0]);
        if(clmparams.clist{ireg}(i) < 0)
            rtmp = fliplr(rtmp);
        end
        iend = istart -1 + nch0*k;
        x(:,istart:iend) = rtmp;
        istart = istart + nch0*k;
    end
    if(ireg>2)
        x(:,nch*k+1) = x(:,1);
    elseif(ireg == 1)
        x(:,1) = [-25;0];
        x(:,nch*k+iextra) = [0;Inf];
        x(:,nch*k+iextra-1) = [25;0];
    else
        x(:,1) = [25;0];
        x(:,nch*k+iextra) = [0;-Inf];
        x(:,nch*k+iextra-1) = [-25;0];
    end
    
    

end