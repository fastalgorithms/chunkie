function uexact = postprocess_uexact_gui(clmparams,targs,targdomain)
    [~,ntarg] = size(targs);
    uexact = zeros(ntarg,1);
    
    if isfield(clmparams,'k')
      k = clmparams.k;
    end
    if isfield(clmparams,'ndomain')
      ndomain = clmparams.ndomain;
    end
    if isfield(clmparams, 'src')
      src = clmparams.src;
    end

    list = cell(1,ndomain);
    for i=1:ndomain
        list{i} = find(targdomain==i);
    end

    
    for i=1:ndomain
        if ~isempty(list{i})    
            j=i+1;
            if j > ndomain
                j = j - ndomain;
            end
            uexact(list{i}) = chnk.helm2d.green(k(i),src(:,j),targs(:,list{i}));
        end
    end

end