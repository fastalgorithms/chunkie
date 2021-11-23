function [x,x1] = get_region_pts_gui_hannah2(chnkr,clmparams,ireg)
    nch = 0;
    k = chnkr(1).k;
    iregstart = 1;
    iregend = length(clmparams.clist{ireg});
    iextra = 1;
    istart = 1;
    idomup = find(clmparams.is_inf == 1);
    idomdown = find(clmparams.is_inf == -1);
    x1  = [];
    
    if(ireg == idomup || ireg == idomdown)
        iregstart = 2;
        iregend = iregend - 1;
        iextra = 3;
        istart = 2;
    end
    if(ireg~=idomdown)
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
        if(ireg~= idomup)
            x(:,nch*k+1) = x(:,1);
        else
            x(:,1) = [-25;0];
            x(:,nch*k+iextra) = [0;Inf];
            x(:,nch*k+iextra-1) = [25;0];
        end
    else
        
        ireg0start = 2;
        ireg0end = 6;
        
        istart = 2;
        iextra = 3;
        for i=ireg0start:ireg0end
            nch = nch + chnkr(abs(clmparams.clist{ireg}(i))).nch;
        end
  
        x = zeros(2,nch*k+iextra);

        for i=ireg0start:ireg0end
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
        x(:,1) = [25;0];
        x(:,nch*k+iextra) = [0;-Inf];
        x(:,nch*k+iextra-1) = [-25;0];
        
        
        
        ireg1start = 7;
        ireg1end = iregend;
        istart = 1;
        nch = 0;
        iextra = 1;
        for i=ireg1start:ireg1end
            nch = nch + chnkr(abs(clmparams.clist{ireg}(i))).nch;
        end
        x1 = zeros(2,nch*k+iextra);
        for i=ireg1start:ireg1end
            rtmp = chnkr(abs(clmparams.clist{ireg}(i))).r;
            nch0 = chnkr(abs(clmparams.clist{ireg}(i))).nch;
            rtmp = reshape(rtmp,[2,k*nch0]);
            if(clmparams.clist{ireg}(i) < 0)
                rtmp = fliplr(rtmp);
            end
            iend = istart -1 + nch0*k;
            x1(:,istart:iend) = rtmp;
            istart = istart + nch0*k;
        end
        
        x1(:,nch*k+1) = x1(:,1);
        
    end
    
    

end