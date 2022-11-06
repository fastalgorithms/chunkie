function [uexact,varargout] = postprocess_uexact_gui(clmparams,targs,targdomain)
    [~,ntarg] = size(targs);
    uexact = zeros(ntarg,1);
    graduexact = zeros(2,ntarg);
    
    if isfield(clmparams,'k')
      k = clmparams.k;
    end
    if isfield(clmparams,'ndomain')
      ndomain = clmparams.ndomain;
    end
    if isfield(clmparams, 'src')
      src = clmparams.src;
    end
     if isfield(clmparams, 'src_in')
      src = clmparams.src_in;
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
            [uexact(list{i}),gtmp] = chnk.helm2d.green(k(i),src(:,j),targs(:,list{i}));
            graduexact(1,list{i}) = reshape(gtmp(:,:,1),[1,length(list{i})]);
            graduexact(2,list{i}) = reshape(gtmp(:,:,2),[1,length(list{i})]);
            
        end
    end
    varargout{1} = graduexact;
end