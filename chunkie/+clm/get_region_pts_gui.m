function [x] = get_region_pts_gui(chnkr,clmparams,ireg)
    nch = 0;
    k = chnkr(1).k;
    for i=1:length(clmparams.clist{ireg})
        nch = nch + chnkr(abs(clmparams.clist{ireg}(i))).nch;
    end
    x = zeros(2,nch*k+1);
    istart = 1;
    for i=1:length(clmparams.clist{ireg})
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
    x(:,nch*k+1) = x(:,1);

end