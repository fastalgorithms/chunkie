function [obj] = balance(obj)
    if (isfield(obj,'vstruc'))
        vstruc = obj.vstruc;
    else
        vstruc = procverts(obj);
    end
    
    echnks = obj.echnks;
    
    nverts = size(obj.verts,2);
    
    for iii=1:nverts
      
    vedge = vstruc{iii}{1};
    vsign = vstruc{iii}{2};
    
    ifdone = false;
    while (~ifdone)
    %%%%%%% find the arclengths of each panel incident to the edge
    
    parcl = zeros([numel(vedge),1]);
    pinds = zeros([numel(vedge),1]);
    
    [xleg16,wleg16,~,~] = lege.exps(16);
    
    for ii=1:numel(vedge)
        if (sign(vsign(ii)) == -1)
            ds =  obj.echnks(vedge(ii)).d(:,:,1);
            pinds(ii) = 1;
        else
            ds =  obj.echnks(vedge(ii)).d(:,:,end);
            pinds(ii) = size(obj.echnks(vedge(ii)).r,3);
        end    
        
        h  =  obj.echnks(vedge(ii)).h(pinds(ii))
        k  =  obj.echnks(vedge(ii)).k;
        if (k ~=16)
             [xleg,wleg,~,~] = lege.exps(k);
        else
            wleg = wleg16;
        end    

        arc = sum(sqrt(ds(1,:).^2+ds(2,:).^2).*wleg'*h);
        parcl(ii) = arc;
          
    end    
    
    amin = min(parcl);
    [amax,ind] = max(parcl);
    ind = ind(1);
    n2ref = floor(log(amax/amin)/log(2));
    
    if (n2ref >0)
        chnkr = obj.echnks(vedge(ind));
        opts = [];
        opts.lvlrfac = 'a';
        for ilev = 1:n2ref
            opts.splitchunks=([pinds(ind)]);
            if (pinds(ind)>1)
                pinds(ind) = pinds(ind)+1;
            end
            chnkr = refine(chnkr,opts);
        end
        obj.echnks(vedge(ind)) = chnkr;
        
    else
        ifdone = true;
    end
    
    end
    end
    
end
