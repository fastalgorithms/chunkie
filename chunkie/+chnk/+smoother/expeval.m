function [cvals] = expeval(czs,ifexps)

    sz = size(czs);
    czs = czs(:);
    cvals= zeros(size(czs));

    ifloc = find(abs(czs) <= 2);
    iffar = find(abs(czs) >2);
    
    if (numel(ifloc))
    cloc = expeval_tayl(czs(ifloc));
    cvals(ifloc) = cloc;
    end
    if (numel(iffar)>0)
    cfar = expeval_lentz(czs(iffar));
    cvals(iffar) = cfar;
    end

    if ((nargin>1) && ifexps)
        cvals = cvals.*exp(-czs);
    end

    cvals = reshape(cvals,sz);

end


function [cvals] = expeval_tayl(czs)

        gamma = 0.577215664901532860606512090082;

        cvals = -gamma-log(czs);
        cpow = czs;
        rfac = 1;

        eps = 1E-24;

        for ii=1:31
            cvals = cvals+cpow.*rfac/ii;
            cpow = -czs.*cpow;
            rfac = rfac/(ii+1);
            % if (max(abs(cpow.*rfac)) < eps)
            %     break;
            % end
        end

         cvals = cvals.*exp(czs);

end


function [cvals] = expeval_lentz(czs)

        nterms = 40;
        cf = 1./czs;
        cd = cf;
        cv = 1/1E-30;

        for ii=1:80
        can = ii;
        cd  = 1./(1+ii*cd);
        cv  = 1 + ii./cv;
        del = cd.*cv;
        cf  = cf.*del;
        cd  = 1./(czs+ii*cd);
        cv  = czs + ii./cv;
        del = cd.*cv;
        cf  = cf.*del;
        % if (max(abs(del-1)) < 1E-14)
        %     break;
        % end
        end

        cvals = cf;

end
