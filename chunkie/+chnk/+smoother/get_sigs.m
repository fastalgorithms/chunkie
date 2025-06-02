function [val,grad] = get_sigs(umesh,rt,sig0,dlam)

        rc = umesh.centroids;
        dlens = umesh.lengths;

        [~,ns] = size(rc);
        [~,nt] = size(rt);

        xc = repmat(rc(1,:),nt,1);
        yc = repmat(rc(2,:),nt,1);
        dl = repmat(dlens,nt,1);

        xt = repmat(rt(1,:).',1,ns);
        yt = repmat(rt(2,:).',1,ns);

        dx = xt-xc;
        dy = yt-yc;
        dr = sqrt(dx.^2+dy.^2);
        dmin = min(dr,[],2);

        dsigj = dl/dlam;
        dexpn = exp(-(dr.^2-dmin.^2)./(2*sig0^2));
        dsignum = sum(dsigj.*dexpn,2);
        dsigden = sum(dexpn,2);
        dsignumx = -sum(dsigj.*dexpn.*(dx./sig0^2),2);
        dsignumy = -sum(dsigj.*dexpn.*(dy./sig0^2),2);
        dsig = dsignum./dsigden;
        dsigx = dsignumx./dsigden;
        dsigy = dsignumy./dsigden;

        val = dsig;

        sz = size(dsig);
        grad = zeros([sz,2]);
        grad(:,:,1) = dsigx;
        grad(:,:,2) = dsigy;
        
        grad = squeeze(grad);

end