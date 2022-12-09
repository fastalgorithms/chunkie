function [ichnk_inds,edgeinds,dens] = proc_graph_regions(cgrph)

    ichnk_inds = [];
    edgeinds = {};
    i0 = 1;
    for ii=1:numel(cgrph.echnks)
        ichnk_inds(ii) = i0;
        chnk_sz = size(cgrph.echnks(ii).r(:,:),2);
        i0 = i0 + chnk_sz;
        edgeinds{ii} = ichnk_inds(ii)-1+(1:chnk_sz);
    end
    ichnk_inds(numel(cgrph.echnks)+1) = i0;

    dens = zeros(ichnk_inds(end)-1,1);
    
    nreg = 0;
    for ireg = 1:numel(cgrph.regions)
        regions = cgrph.regions{ireg};
        for icomp = 1:numel(regions)
            component = regions{icomp};
            nreg = nreg + 1;
           	for icurve = 1:numel(component)
            	enums = component{icurve};
                for ichnk = 1:numel(enums)
                    inds = edgeinds{abs(enums(ichnk))};
                    
                    dens(inds) = dens(inds) + sign(enums(ichnk))*2^(nreg-1);
                end
            end
        end
    end                
end

